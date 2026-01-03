#ifndef REMESH_BOTSCH
#define REMESH_BOTSCH

#include <Eigen/Core>
#include <Eigen/src/Core/Matrix.h>

void remesh_botsch(const Eigen::MatrixXd &Vattrs_in,
                   const Eigen::MatrixX3i &F_in,
                   const Eigen::VectorXi &Vselection_in,
                   const Eigen::VectorXd &Vtargetlen_in,
                   double selection_threshold, int iters, bool smooth,
                   bool project, int verbose, Eigen::MatrixXd &Vattrs_out,
                   Eigen::MatrixX3i &F_out, Eigen::VectorXi &Vselection_out,
                   Eigen::VectorXi &Fi_containing_V_proj_out,
                   Eigen::VectorXi &new2oldFi_out);

// overload with constant targetlen
void remesh_botsch(const Eigen::MatrixXd &Vattrs_in,
                   const Eigen::MatrixX3i &F_in,
                   const Eigen::VectorXi &Vselection_in, double targetlen,
                   double selection_threshold, int iters, bool smooth,
                   bool project, int verbose, Eigen::MatrixXd &Vattrs_out,
                   Eigen::MatrixX3i &F_out, Eigen::VectorXi &Vselection_out,
                   Eigen::VectorXi &Fi_containing_V_proj_out,
                   Eigen::VectorXi &new2oldFi_out);

void remesh_botsch_adaptive(
    const Eigen::MatrixXd &Vattrs_in, const Eigen::MatrixX3i &F_in,
    const Eigen::VectorXi &Vselection_in, double epsilon, bool adaptive,
    double selection_threshold, int iters, bool smooth, bool project,
    int verbose, Eigen::MatrixXd &Vattrs_out, Eigen::MatrixX3i &F_out,
    Eigen::VectorXi &Vselection_out, Eigen::VectorXi &Fi_containing_V_proj_out,
    Eigen::VectorXi &new2oldFi_out);

// helper for finding barycentric hits with Fi_containing_V_proj
void barycentric_interp_on_Fi_containing_V_proj(
    const Eigen::MatrixXd &Vattrs_orig_in, const Eigen::MatrixX3i &F_orig_in,
    const Eigen::MatrixX3d &V_proj_in,
    const Eigen::VectorXi &Fi_containing_V_proj_in,
    Eigen::MatrixX3d &barycentric_coordinates_out,
    Eigen::MatrixXd &Vattrs_barylerped_out);

#endif
