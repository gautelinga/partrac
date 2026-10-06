#include "SimplexFreqInterpol.hpp"
#include <iostream>

// Out of the evaluation's unit, which it would otherwise grow

template<typename Cell>
void SimplexFreqInterpol<Cell>::freeze(const double t)
{
  update(t);
  const FreqWeights wf = fill_weights(t);
  const auto sum = [&](std::vector<std::vector<double>>& nodes){
    if (nodes.empty()) return;
    std::vector<double> s(nodes[0].size(), 0.);
    for (std::size_t k = 0; k < nodes.size(); ++k)
      for (std::size_t i = 0; i < s.size(); ++i) s[i] += wf.w[k]*nodes[k][i];
    nodes.assign(1, std::move(s));
  };
  sum(u_nodes_);
  sum(p_nodes_);
  fs.make_steady();
  // Weights kept for the components before
  new_id();
  std::cout << "Fields frozen at t = " << t << std::endl;
}

template void SimplexFreqInterpol<Triangle>::freeze(const double t);
template void SimplexFreqInterpol<Tet>::freeze(const double t);
