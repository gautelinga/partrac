// The held evaluation, in a unit of its own: apart from the located evaluation, whose inlining it would change
#include "SplitInterpol_eval.hpp"

template void SplitInterpol<Triangle>::hold(const SubRegion&, Held&) const;
template void SplitInterpol<Tet>::hold(const SubRegion&, Held&) const;
template void SplitInterpol<Triangle>::held_motion(int, const std::array<double, 4>&, double, const Held&, PointValues&);
template void SplitInterpol<Tet>::held_motion(int, const std::array<double, 4>&, double, const Held&, PointValues&);
template void SplitInterpol<Triangle>::held_velocity(int, const std::array<double, 4>&, double, const Held&, PointValues&);
template void SplitInterpol<Tet>::held_velocity(int, const std::array<double, 4>&, double, const Held&, PointValues&);
