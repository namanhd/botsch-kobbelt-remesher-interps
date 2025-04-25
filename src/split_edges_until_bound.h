#ifndef SPLIT_EDGES_UNTIL_BOUND
#define SPLIT_EDGES_UNTIL_BOUND

#include <Eigen/Core>

// old signature without templating
// void split_edges_until_bound(Eigen::MatrixXd &V_etc, Eigen::MatrixXi &F,
// double selthresh, Eigen::VectorXi &feature,
// Eigen::VectorXd &high, Eigen::VectorXd &low);

template <typename DerivedV_etc, typename DerivedF, typename DerivedFeature,
          typename DerivedHigh, typename DerivedLow>
void split_edges_until_bound(Eigen::PlainObjectBase<DerivedV_etc> &V_etc,
                             Eigen::PlainObjectBase<DerivedF> &F,
                             double selthresh,
                             Eigen::PlainObjectBase<DerivedFeature> &feature,
                             Eigen::PlainObjectBase<DerivedHigh> &high,
                             Eigen::PlainObjectBase<DerivedLow> &low);

// header-only needed to write eigen functions in this templated manner
// to prevent copying and to have the most generic, least-copy eigen code
#include "split_edges_until_bound.cpp"
#endif
