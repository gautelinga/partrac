#ifndef __PERF_WINDOW_HPP
#define __PERF_WINDOW_HPP

// Counters of a surrounding `perf stat --delay=-1 --control=fifo:ctl,ack` run,
// on around the time loop only: PARTRAC_PERF_CTL and PARTRAC_PERF_ACK name the
// two fifos; unset, nothing happens

#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <unistd.h>

class PerfWindow {
public:
  PerfWindow(){
    const char* ctl = std::getenv("PARTRAC_PERF_CTL");
    const char* ack = std::getenv("PARTRAC_PERF_ACK");
    if (!ctl || !ack) return;
    // Read-write, so neither open waits for the other end
    ctl_ = open(ctl, O_RDWR);
    ack_ = open(ack, O_RDWR);
  }
  ~PerfWindow(){
    if (ctl_ >= 0) close(ctl_);
    if (ack_ >= 0) close(ack_);
  }
  void enable() { send("enable\n"); }
  void disable() { send("disable\n"); }
private:
  void send(const char* cmd){
    if (ctl_ < 0 || ack_ < 0) return;
    if (write(ctl_, cmd, std::strlen(cmd)) < 0) return;
    char buf[16];
    if (read(ack_, buf, sizeof(buf)) < 0) return;   // perf answers "ack"
  }
  int ctl_ = -1;
  int ack_ = -1;
};

#endif
