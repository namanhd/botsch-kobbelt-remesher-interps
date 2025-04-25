#ifndef COLLAPSE_EDGES
#define COLLAPSE_EDGES

#include <Eigen/Core>

void collapse_edges(Eigen::MatrixXd &V_etc, Eigen::MatrixXi &F,
                    double selthresh, Eigen::VectorXi &feature,
                    Eigen::VectorXd &high, Eigen::VectorXd &low);

#endif
