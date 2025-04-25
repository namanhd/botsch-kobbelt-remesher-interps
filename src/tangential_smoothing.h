#ifndef TANGENTIAL_RELAXATION
#define TANGENTIAL_RELAXATION

#include <Eigen/Core>

// old signature without templating
// void tangential_smoothing(Eigen::MatrixXd &V_etc, Eigen::MatrixXi &F,
// double selthresh, Eigen::VectorXi &feature,
// Eigen::VectorXd &lambda,
// const Eigen::VectorXd &sizingField,
// bool use_sizing_field_weighted_avg);

template <typename DerivedV_etc, typename DerivedF, typename DerivedFeature,
          typename DerivedLambda, typename DerivedSizingField>
void tangential_smoothing(
    Eigen::MatrixBase<DerivedV_etc> &V_etc, Eigen::MatrixBase<DerivedF> &F,
    double selthresh, Eigen::MatrixBase<DerivedFeature> &feature,
    Eigen::MatrixBase<DerivedLambda> &lambda,
    const Eigen::MatrixBase<DerivedSizingField> &sizingField,
    bool use_sizing_field_weighted_avg);

// header-only needed to write eigen functions in this templated manner
// to prevent copying and to have the most generic, least-copy eigen code
#include "tangential_smoothing.cpp"
#endif
