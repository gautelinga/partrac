#ifndef __STATS_COLUMNS_HPP
#define __STATS_COLUMNS_HPP

#include <fstream>
#include <vector>

// A statistic is a name and a value together, so the header and the rows are
// two readings of one list
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
