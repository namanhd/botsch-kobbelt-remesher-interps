#ifndef REMESH_BOTSCH
#define REMESH_BOTSCH

#include <Eigen/Core>

void remesh_botsch(const Eigen::MatrixXd &Vattrs_in,
                   const Eigen::MatrixXi &F_in,
                   const Eigen::VectorXi &Vselection_in,
                   const Eigen::VectorXd &Vtargetlen_in,
                   double selection_threshold, int iters, bool project,
                   int verbose, Eigen::MatrixXd &Vattrs_out,
                   Eigen::MatrixXi &F_out, Eigen::VectorXi &Vselection_out);

// overload with constant targetlen
void remesh_botsch(const Eigen::MatrixXd &Vattrs_in,
                   const Eigen::MatrixXi &F_in,
                   const Eigen::VectorXi &Vselection_in, double targetlen,
                   double selection_threshold, int iters, bool project,
                   int verbose, Eigen::MatrixXd &Vattrs_out,
                   Eigen::MatrixXi &F_out, Eigen::VectorXi &Vselection_out);

void remesh_botsch_adaptive(const Eigen::MatrixXd &Vattrs_in,
                            const Eigen::MatrixXi &F_in,
                            const Eigen::VectorXi &Vselection_in,
                            double epsilon, bool adaptive,
                            double selection_threshold, int iters, bool project,
                            int verbose, Eigen::MatrixXd &Vattrs_out,
                            Eigen::MatrixXi &F_out,
                            Eigen::VectorXi &Vselection_out);

#endif
