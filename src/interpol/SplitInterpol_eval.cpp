#include "SplitInterpol_eval.hpp"

template void SplitInterpol<Triangle>::evaluate(const Vector3d&, double, const CellPos&, PointValues&);
template void SplitInterpol<Triangle>::evaluate_motion(const Vector3d&, double, const CellPos&, PointValues&);
template void SplitInterpol<Tet>::evaluate(const Vector3d&, double, const CellPos&, PointValues&);
template void SplitInterpol<Tet>::evaluate_motion(const Vector3d&, double, const CellPos&, PointValues&);
