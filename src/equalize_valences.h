#ifndef EQUALIZE_VALENCES
#define EQUALIZE_VALENCES

#include <Eigen/Core>

void equalize_valences(Eigen::MatrixXd &V_etc, Eigen::MatrixXi &F,
                       double selthresh, Eigen::VectorXi &feature);

#endif
