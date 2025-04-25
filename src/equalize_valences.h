#ifndef EQUALIZE_VALENCES
#define EQUALIZE_VALENCES

#include <Eigen/Core>

// old signature without templating
// void equalize_valences(Eigen::MatrixXd &V_etc, Eigen::MatrixXi &F,
// double selthresh, Eigen::VectorXi &feature);

template <typename DerivedV_etc, typename DerivedF, typename DerivedFeature>
void equalize_valences(Eigen::PlainObjectBase<DerivedV_etc> &V_etc,
                       Eigen::PlainObjectBase<DerivedF> &F, double selthresh,
                       Eigen::PlainObjectBase<DerivedFeature> &feature);

// header-only needed to write eigen functions in this templated manner
// to prevent copying and to have the most generic, least-copy eigen code
#include "equalize_valences.cpp"
#endif
