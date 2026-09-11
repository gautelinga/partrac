#ifndef __STATS_COLUMNS_HPP
#define __STATS_COLUMNS_HPP

#include <fstream>
#include <vector>

// A statistic is a name and a value together, so a header and its rows are two
// readings of one list rather than two functions kept in step by hand. They had
// drifted five times before this.
//
// What is shared here is the mechanism, not the columns: partrac reports a
// mesh, the experimental apps report a particle set, and the tensor tracers
// report their own thing. Those are three different statistics and merging them
// would be a modelling decision, not a refactor. Drift is what they had in
// common, and this is what removes it.
struct StatsColumn {
  const char* name;
  double value;
  bool is_count = false;   // written as an integer, not in the double format
};

inline void write_stats_header(std::ofstream &statfile,
                               const std::vector<StatsColumn>& cols){
  statfile << "# ";
  for (const auto& col : cols)
    statfile << col.name << "\t";
  statfile << std::endl;
}

inline void write_stats_row(std::ofstream &statfile,
                            const std::vector<StatsColumn>& cols){
  for (const auto& col : cols){
    if (col.is_count)
      statfile << static_cast<unsigned long long>(col.value) << "\t";
    else
      statfile << col.value << "\t";
  }
  statfile << std::endl;
}

#endif
