// The stamp buffer both mesh loaders blend a field through: which stamps it
// reads, and which it must not. A stamp is a file read, seconds of it on a
// large mesh, and a run that steps forward asks for the same bracket over and
// over and then for one shifted by a stamp, so a buffer that re-reads what it
// already holds is silent except in the wall clock. What travels with a stamp
// matters as much: the near-wall tolerance of the XDMF loader is computed from
// a stamp's velocity, and a swap that left it behind would give a wall the
// wrong rule. The reads are counted through the callable, so a re-read shows.
#include <catch2/catch.hpp>

#include <string>
#include <vector>

#include "stamp_buffer.hpp"

namespace {

// A stamp: one value and something computed from it, as a loader's near-wall
// tolerance is
struct Stamp {
  double value = 0.;
  double note = 0.;
};

// Reads a stamp named by an int, counting what it was asked to read
struct Reader {
  std::vector<int> read;
  void operator()(const int key, Stamp& s){
    read.push_back(key);
    s.value = 10.*key;
    s.note = 100.*key;
  }
  std::size_t count() const { return read.size(); }
};

using Buffer = partrac::StampBuffer<Stamp, int>;

}  // namespace

TEST_CASE("The stamp buffer reads a stamp once", "[stamp_buffer]") {
  Buffer buf;
  Reader r;

  // The first bracket: both stamps are new
  auto fill = buf.load(0, 1, r);
  REQUIRE(r.read == std::vector<int>{0, 1});
  REQUIRE(fill.first == partrac::StampFill::Read);
  REQUIRE(fill.second == partrac::StampFill::Read);
  REQUIRE(buf.prev().value == 0.);
  REQUIRE(buf.next().value == 10.);

  // The same bracket again, as every step within it asks for it
  fill = buf.load(0, 1, r);
  REQUIRE(r.count() == 2);
  REQUIRE(fill.first == partrac::StampFill::Held);
  REQUIRE(fill.second == partrac::StampFill::Held);
  REQUIRE(buf.prev().value == 0.);
  REQUIRE(buf.next().value == 10.);

  // One stamp on: what was next is now previous, and only the new one is read
  fill = buf.load(1, 2, r);
  REQUIRE(r.read == std::vector<int>{0, 1, 2});
  REQUIRE(fill.first == partrac::StampFill::Swapped);
  REQUIRE(fill.second == partrac::StampFill::Read);
  REQUIRE(buf.prev().value == 10.);
  REQUIRE(buf.next().value == 20.);

  // A bracket neither buffer holds costs both stamps
  fill = buf.load(7, 8, r);
  REQUIRE(r.read == std::vector<int>{0, 1, 2, 7, 8});
  REQUIRE(fill.first == partrac::StampFill::Read);
  REQUIRE(fill.second == partrac::StampFill::Read);
}

TEST_CASE("A step back reads only the stamp it has not got", "[stamp_buffer]") {
  // Stepping back a stamp, the previous stamp becomes the next one: it moves
  // to the next stamp's buffer instead of being read again
  Buffer buf;
  Reader r;
  buf.load(1, 2, r);
  const auto fill = buf.load(0, 1, r);
  REQUIRE(r.read == std::vector<int>{1, 2, 0});
  REQUIRE(fill.first == partrac::StampFill::Read);
  REQUIRE(fill.second == partrac::StampFill::Held);
  REQUIRE(buf.prev().value == 0.);
  REQUIRE(buf.next().value == 10.);
  REQUIRE(buf.next().note == 100.);

  // and a bracket given the other way round is one swap, no read
  buf.load(1, 0, r);
  REQUIRE(r.count() == 3);
  REQUIRE(buf.prev().value == 10.);
  REQUIRE(buf.next().value == 0.);
}

TEST_CASE("A stamp's own values travel with it", "[stamp_buffer]") {
  Buffer buf;
  Reader r;
  buf.load(3, 4, r);
  const Stamp* const held_4 = &buf.next();
  REQUIRE(buf.next().note == 400.);

  // The swap moves the stamp, not a copy of it, and what was computed from it
  // goes along
  const auto fill = buf.load(4, 5, r);
  REQUIRE(fill.first == partrac::StampFill::Swapped);
  REQUIRE(buf.prev().note == 400.);
  REQUIRE(buf.next().note == 500.);
  REQUIRE(&buf.prev() != held_4);   // the buffers, not their contents, are fixed
  REQUIRE(r.read == std::vector<int>{3, 4, 5});
}

TEST_CASE("A single stamp is aliased, not copied", "[stamp_buffer]") {
  Buffer buf;
  Reader r;

  // A field with one stamp brackets it with itself: one read, one buffer
  const auto fill = buf.load(2, 2, r);
  REQUIRE(r.read == std::vector<int>{2});
  REQUIRE(fill.first == partrac::StampFill::Read);
  REQUIRE(fill.second == partrac::StampFill::Aliased);
  REQUIRE(buf.aliased());
  REQUIRE(&buf.prev() == &buf.next());
  REQUIRE(buf.next().value == 20.);

  // and a bracket of two stamps after it is no longer aliased
  buf.load(2, 3, r);
  REQUIRE_FALSE(buf.aliased());
  REQUIRE(r.read == std::vector<int>{2, 3});
}

TEST_CASE("The first stamp may be held by a loader's setup", "[stamp_buffer]") {
  Buffer buf;
  Reader r;

  // A loader that builds its tables from the first stamp keeps its values
  buf.a().value = 10.;
  buf.a().note = 100.;
  buf.hold_a(1);
  REQUIRE(buf.key_a() == 1);

  const auto fill = buf.load(1, 2, r);
  REQUIRE(r.read == std::vector<int>{2});
  REQUIRE(fill.first == partrac::StampFill::Held);
  REQUIRE(buf.prev().note == 100.);
}

TEST_CASE("An empty buffer is not the stamp its key defaults to", "[stamp_buffer]") {
  // Nothing is held before the first load, and a key that is zero, or empty,
  // is a stamp like any other
  Buffer buf;
  Reader r;
  buf.load(0, 0, r);
  REQUIRE(r.read == std::vector<int>{0});

  partrac::StampBuffer<Stamp, std::string> named;
  std::vector<std::string> read;
  auto by_name = [&](const std::string& key, Stamp& s){ read.push_back(key); s.value = 1.; };
  named.load("", "", by_name);
  REQUIRE(read == std::vector<std::string>{""});
}

TEST_CASE("The bracket is loaded once, and kept past t_max", "[stamp_buffer]") {
  const partrac::StampTimes first{0., 1.};
  const partrac::StampTimes second{1., 2.};

  // Nothing held: the first update loads whatever the time asks for, even at
  // or past the end of the field data
  REQUIRE(partrac::stamp_reload(false, 0., 2., {0., 0.}, first));
  REQUIRE(partrac::stamp_reload(false, 2., 2., {0., 0.}, second));
  REQUIRE(partrac::stamp_reload(false, 9., 2., {0., 0.}, second));

  // The bracket held is the bracket wanted: nothing to do
  REQUIRE_FALSE(partrac::stamp_reload(true, 0.5, 2., first, first));

  // A step into the next bracket, within the field data: load it
  REQUIRE(partrac::stamp_reload(true, 1.5, 2., first, second));

  // At and past t_max the last bracket is kept, so a run that ends there still
  // has two stamps to take a rate from
  REQUIRE_FALSE(partrac::stamp_reload(true, 2., 2., first, second));
  REQUIRE_FALSE(partrac::stamp_reload(true, 3., 2., first, second));
}
