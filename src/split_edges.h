#ifndef SPLIT_EDGES
#define SPLIT_EDGES

#include <Eigen/Core>
#include <vector>

// old signature without templating
// void split_edges(Eigen::MatrixXd &V_etc, Eigen::MatrixXi &F,
//                  Eigen::MatrixXi &E0, Eigen::MatrixXi &uE,
//                  Eigen::VectorXi &EMAP0, std::vector<std::vector<int>> &uE2E,
//                  Eigen::VectorXd &high, Eigen::VectorXd &low,
//                  const std::vector<int> &edges_to_split);

template <typename DerivedV_etc, typename DerivedF, typename DerivedE0,
          typename DeriveduE, typename DerivedEMAP0, typename uE2EType,
          typename DerivedHigh, typename DerivedLow>
void split_edges(Eigen::PlainObjectBase<DerivedV_etc> &V_etc,
                 Eigen::PlainObjectBase<DerivedF> &F,
                 Eigen::PlainObjectBase<DerivedE0> &E0,
                 Eigen::PlainObjectBase<DeriveduE> &uE,
                 Eigen::PlainObjectBase<DerivedEMAP0> &EMAP0,
                 std::vector<std::vector<uE2EType>> &uE2E,
                 Eigen::PlainObjectBase<DerivedHigh> &high,
                 Eigen::PlainObjectBase<DerivedLow> &low,
                 const std::vector<uE2EType> &edges_to_split);

#include "split_edges.cpp"
#endif
