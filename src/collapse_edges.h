#ifndef COLLAPSE_EDGES
#define COLLAPSE_EDGES

#include <Eigen/Core>

// old signature without templating
// void collapse_edges(Eigen::MatrixXd &V_etc, Eigen::MatrixXi &F,
//   double selthresh, Eigen::VectorXi &feature,
//   Eigen::VectorXd &high, Eigen::VectorXd &low);

template <typename DerivedV_etc, typename DerivedF, typename DerivedFeature,
          typename DerivedHigh, typename DerivedLow, typename Derivednew2oldFi>
void collapse_edges(Eigen::MatrixBase<DerivedV_etc> &V_etc,
                    Eigen::MatrixBase<DerivedF> &F, double selthresh,
                    Eigen::MatrixBase<DerivedFeature> &feature,
                    Eigen::MatrixBase<DerivedHigh> &high,
                    Eigen::MatrixBase<DerivedLow> &low,
                    Eigen::MatrixBase<Derivednew2oldFi> &new2oldFi);

// header-only needed to write eigen functions in this templated manner
// to prevent copying and to have the most generic, least-copy eigen code
#include "collapse_edges.cpp"
#endif
