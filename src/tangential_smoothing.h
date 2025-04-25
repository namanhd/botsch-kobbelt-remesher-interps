#ifndef TANGENTIAL_RELAXATION
#define TANGENTIAL_RELAXATION



#include <Eigen/Core>

void tangential_smoothing(Eigen::MatrixXd &V_etc, Eigen::MatrixXi &F, double selthresh, Eigen::VectorXi &feature, Eigen::VectorXd &lambda, 
    const Eigen::VectorXd &sizingField, bool use_sizing_field_weighted_avg);


#endif
