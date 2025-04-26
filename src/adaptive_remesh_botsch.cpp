#include "adaptive_remesh_botsch.h"
#include "collapse_edges.h"
#include "equalize_valences.h"
#include "igl/barycentric_coordinates.h"
#include "igl/barycentric_interpolation.h"
#include "tangential_smoothing.h"
#include <Eigen/Core>
// #include <igl/is_edge_manifold.h>
// #include <igl/writeOBJ.h>
#include "split_edges_until_bound.h"
// #include <igl/unique_edge_map.h>
// #include <igl/edge_flaps.h>
// #include <igl/circulation.h>
// #include <igl/remove_duplicate_vertices.h>
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Core/util/Constants.h>
#include <igl/avg_edge_length.h>
#include <igl/cotmatrix.h>
#include <igl/gaussian_curvature.h>
#include <igl/invert_diag.h>
#include <igl/massmatrix.h>
#include <igl/point_mesh_squared_distance.h>
#include <iostream>

template <typename DerivedV_etc, typename DerivedF, typename DerivedSizingField>
void calc_adaptive_sizing_field(
    Eigen::MatrixBase<DerivedV_etc> &V_etc, Eigen::MatrixBase<DerivedF> &F,
    Eigen::PlainObjectBase<DerivedSizingField> &sizingField, double epsilon,
    bool adaptive, double my_min_adaptive_sizing,
    double my_max_adaptive_sizing) {
  if (adaptive) {
    Eigen::VectorXd epsilon_vec =
        Eigen::VectorXd::Constant(V_etc.rows(), epsilon);
    Eigen::MatrixXd K_tot; // maximum absolute curvature
    Eigen::VectorXd K;     // gaussian curvature
    // Compute integral of Gaussian curvature
    igl::gaussian_curvature(V_etc.template leftCols<3>(), F, K);
    // Compute mass matrix
    Eigen::SparseMatrix<double> L, M, Minv;
    igl::massmatrix(V_etc.template leftCols<3>(), F,
                    igl::MASSMATRIX_TYPE_VORONOI, M);
    igl::invert_diag(M, Minv);
    // Divide by area to get integral average
    K = (Minv * K).eval();
    // mean curvature
    Eigen::MatrixXd HN, H;
    igl::cotmatrix(V_etc.template leftCols<3>(), F, L);
    HN = -Minv * (L * V_etc.template leftCols<3>());
    H = HN.rowwise().norm().array() / 4; // up to sign.
    Eigen::VectorXd delta = H.array().square() - K.array();
    for (int i = 0; i < delta.size(); i++) // clip delta to 0 to avoid nans
    {
      if (delta[i] < 0) {
        delta(i) = 0;
      }
    }
    K_tot = H.array() + delta.array().sqrt();
    K_tot = K_tot.matrix();
    sizingField = (6 * epsilon_vec.array() * K_tot.array().inverse() -
                   3 * epsilon_vec.array().square())
                      .sqrt(); // this could potentially become nan. Authors did
                               // not discuss this.
    // std::cout << "max: " <<
    // sizingField.array().isNaN().select(0,sizingField).maxCoeff()
    // << "\n"; std::cout << "min: " <<
    // sizingField.array().isNaN().select(0,sizingField).minCoeff() << "\n";
    for (int i = 0; i < sizingField.size();
         i++) // clip sizingField if it is below or above values
    {

      // 	// these are yoterel's values
      // 	// double my_max = 4;
      // 	// double my_min = 0.01;
      if ((sizingField(i) >
           my_max_adaptive_sizing)) // || std::isnan(sizingField(i))
      {
        sizingField(i) = my_max_adaptive_sizing;
      }
      if ((sizingField(i) < my_min_adaptive_sizing)) {
        sizingField(i) = my_min_adaptive_sizing;
      }
      // std::cout << "found nan in sizing field. i: " << i << "\n";
      // std::cout << "ktot: " << K_tot.array()(i) << "\n";
      // std::cout << "epsilon: " << epsilon.array()(i) << "\n";
      // std::cout << "H: " << H.array()(i) << "\n";
      // std::cout << "K: " << K.array()(i) << "\n";
    }
    sizingField = sizingField.matrix();
  } else {
    sizingField = Eigen::VectorXd::Constant(V_etc.rows(), epsilon);
  }
}

void remesh_botsch(const Eigen::MatrixXd &Vattrs_in,
                   const Eigen::MatrixX3i &F_in,
                   const Eigen::VectorXi &Vselection_in,
                   const Eigen::VectorXd &Vtargetlen_in,
                   double selection_threshold, int iters, bool project,
                   int verbose, Eigen::MatrixXd &Vattrs_out,
                   Eigen::MatrixX3i &F_out, Eigen::VectorXi &Vselection_out,
                   Eigen::VectorXi &Fi_containing_V_proj_out) {
  Eigen::MatrixXd V0;
  Eigen::MatrixX3i F0;
  Eigen::VectorXd high, low, lambda, sizingField;
  // high = Vtargetlen_in * 1.4;
  // low = Vtargetlen_in * 0.7;

  F0 = F_in;
  F_out = F0;
  V0 = Vattrs_in.leftCols<3>();

  Eigen::VectorXi feature;
  feature.resize(0);

  // stack Vselection_in and Vtargetlen_in onto Vattrs_in, forming V_etc that
  // will be interpolated throughout
  int n_attribs = Vattrs_in.cols();
  Eigen::MatrixXd V_etc(Vattrs_in.rows(), n_attribs + 2);
  V_etc << Vattrs_in, Vtargetlen_in, Vselection_in.cast<double>();
  if (verbose) {
    std::cout << "V_etc shape " << V_etc.rows() << " " << V_etc.cols() << "\n";
  }

  // Iterate the four steps
  for (int i = 0; i < iters; i++) {
    if (verbose) {
      std::cout << "iter: " << i << "/" << iters << "\n";
    }
    // grab the latest interpd sizingField from the second-to-last V_etc coord
    sizingField = V_etc(Eigen::all, Eigen::last - Eigen::fix<1>);
    high = 1.4 * sizingField;
    low = 0.7 * sizingField;

    if (verbose) {
      std::cout << "splitting...\n";
    }
    split_edges_until_bound(V_etc, F_out, selection_threshold, feature, high,
                            low); // Split

    if (verbose) {
      std::cout << "collapsing...\n";
    }
    collapse_edges(V_etc, F_out, selection_threshold, feature, high,
                   low); // Collapse

    int new_n_verts = V_etc.rows();

    // refresh sizingField for latest V_etc size, for  tangential_smoothing
    sizingField = V_etc(Eigen::all, Eigen::last - Eigen::fix<1>);

    if (verbose) {
      std::cout << "flipping...\n";
    }
    equalize_valences(V_etc, F_out, selection_threshold, feature); // Flip

    lambda = Eigen::VectorXd::Constant(new_n_verts, 1.0);
    if (verbose) {
      std::cout << "smoothing...\n";
    }
    tangential_smoothing(V_etc, F_out, selection_threshold, feature, lambda,
                         sizingField, true); // Smooth

    if (project) {
      if (verbose) {
        std::cout << "projecting..." << "\n";
      }
      Eigen::VectorXd sqrD;
      Eigen::MatrixX3d V_projected;
      igl::point_mesh_squared_distance(V_etc.leftCols(3).eval(), V0, F0, sqrD,
                                       Fi_containing_V_proj_out,
                                       V_projected); // Project
      V_etc.leftCols<3>() = V_projected;
    }

    if (verbose) {
      std::cout << "finished iter, faces amount: " << F_out.rows() << "\n";
    }
  }
  // return in Vattrs_out, removing the appended Vselection and Vtargetlen in
  // V_etc
  Vattrs_out = V_etc(Eigen::all, Eigen::seq(0, Eigen::last - Eigen::fix<2>));
  // return in Vselection_out
  Vselection_out = (V_etc.rightCols<1>().array() >= selection_threshold)
                       .cast<int>()
                       .matrix();
}

// overload with constant targetlen, does not need to append the targetlen field
// onto V_etc
void remesh_botsch(const Eigen::MatrixXd &Vattrs_in,
                   const Eigen::MatrixX3i &F_in,
                   const Eigen::VectorXi &Vselection_in, double targetlen,
                   double selection_threshold, int iters, bool project,
                   int verbose, Eigen::MatrixXd &Vattrs_out,
                   Eigen::MatrixX3i &F_out, Eigen::VectorXi &Vselection_out,
                   Eigen::VectorXi &Fi_containing_V_proj_out) {
  Eigen::MatrixXd V0;
  Eigen::MatrixX3i F0;
  Eigen::VectorXd high, low, lambda, sizingField;

  F0 = F_in;
  F_out = F0;
  V0 = Vattrs_in.leftCols<3>();

  Eigen::VectorXi feature;
  feature.resize(0);

  // stack Vselection_in onto Vattrs_in, forming V_etc that will be interpolated
  // throughout
  int n_attribs = Vattrs_in.cols();
  Eigen::MatrixXd V_etc(Vattrs_in.rows(), n_attribs + 1);
  V_etc << Vattrs_in, Vselection_in.cast<double>();
  if (verbose) {
    std::cout << "V_etc shape " << V_etc.rows() << " " << V_etc.cols() << "\n";
  }

  // Iterate the four steps
  for (int i = 0; i < iters; i++) {
    if (verbose) {
      std::cout << "iter: " << i << "/" << iters << "\n";
    }

    sizingField = Eigen::VectorXd::Constant(V_etc.rows(), targetlen);
    high = 1.4 * sizingField;
    low = 0.7 * sizingField;

    if (verbose) {
      std::cout << "splitting...\n";
    }
    split_edges_until_bound(V_etc, F_out, selection_threshold, feature, high,
                            low); // Split

    if (verbose) {
      std::cout << "collapsing...\n";
    }
    collapse_edges(V_etc, F_out, selection_threshold, feature, high,
                   low); // Collapse

    int new_n_verts = V_etc.rows();
    // update this for the latest V_etc size which won't change anymore for the
    // rest of this loop iteration
    sizingField = Eigen::VectorXd::Constant(new_n_verts, targetlen);

    if (verbose) {
      std::cout << "flipping...\n";
    }
    equalize_valences(V_etc, F_out, selection_threshold, feature); // Flip

    lambda = Eigen::VectorXd::Constant(new_n_verts, 1.0);
    if (verbose) {
      std::cout << "smoothing...\n";
    }
    tangential_smoothing(V_etc, F_out, selection_threshold, feature, lambda,
                         sizingField,
                         true); // Smooth

    if (project) {
      if (verbose) {
        std::cout << "projecting..." << "\n";
      }
      Eigen::VectorXd sqrD;
      Eigen::MatrixX3d V_projected;
      igl::point_mesh_squared_distance(V_etc.leftCols(3).eval(), V0, F0, sqrD,
                                       Fi_containing_V_proj_out,
                                       V_projected); // Project
      V_etc.leftCols<3>() = V_projected;
    }

    if (verbose) {
      std::cout << "finished iter, faces amount: " << F_out.rows() << "\n";
    }
  }
  // return in Vattrs_out, removing the appended Vselection in V_etc
  Vattrs_out = V_etc(Eigen::all, Eigen::seq(0, Eigen::last - Eigen::fix<1>));
  // return in Vselection_out
  Vselection_out = (V_etc.rightCols<1>().array() >= selection_threshold)
                       .cast<int>()
                       .matrix();
}

// this really should be an overload...
void remesh_botsch_adaptive(const Eigen::MatrixXd &Vattrs_in,
                            const Eigen::MatrixX3i &F_in,
                            const Eigen::VectorXi &Vselection_in,
                            double epsilon, bool adaptive,
                            double selection_threshold, int iters, bool project,
                            int verbose, Eigen::MatrixXd &Vattrs_out,
                            Eigen::MatrixX3i &F_out,
                            Eigen::VectorXi &Vselection_out,
                            Eigen::VectorXi &Fi_containing_V_proj_out) {
  Eigen::MatrixXd V0;
  Eigen::MatrixX3i F0;
  Eigen::VectorXd high, low, lambda, sizingField;

  F0 = F_in;
  F_out = F0;
  V0 = Vattrs_in.leftCols<3>();

  Eigen::VectorXi feature;
  feature.resize(0);

  // stack Vselection_in onto Vattrs_in, forming V_etc that will be interpolated
  // throughout
  int n_attribs = Vattrs_in.cols();
  Eigen::MatrixXd V_etc(Vattrs_in.rows(), n_attribs + 1);
  V_etc << Vattrs_in, Vselection_in.cast<double>();
  if (verbose) {
    std::cout << "V_etc shape " << V_etc.rows() << " " << V_etc.cols() << "\n";
  }

  // Iterate the four steps
  for (int i = 0; i < iters; i++) {
    if (verbose) {
      std::cout << "iter: " << i << "/" << iters << "\n";
    }

    calc_adaptive_sizing_field(V_etc, F_out, sizingField, epsilon, adaptive,
                               0.01, 4.0);
    high = 1.4 * sizingField;
    low = 0.7 * sizingField;

    if (verbose) {
      std::cout << "splitting...\n";
    }
    split_edges_until_bound(V_etc, F_out, selection_threshold, feature, high,
                            low); // Split

    if (verbose) {
      std::cout << "collapsing...\n";
    }
    collapse_edges(V_etc, F_out, selection_threshold, feature, high,
                   low); // Collapse

    // update this for the latest V_etc size which won't change anymore for the
    // rest of this loop iteration
    calc_adaptive_sizing_field(V_etc, F_out, sizingField, epsilon, adaptive,
                               0.01, 4.0);

    int new_n_verts = V_etc.rows();

    if (verbose) {
      std::cout << "flipping...\n";
    }
    equalize_valences(V_etc, F_out, selection_threshold, feature); // Flip

    lambda = Eigen::VectorXd::Constant(new_n_verts, 1.0);
    if (verbose) {
      std::cout << "smoothing...\n";
    }
    tangential_smoothing(V_etc, F_out, selection_threshold, feature, lambda,
                         sizingField,
                         true); // Smooth

    if (project) {
      if (verbose) {
        std::cout << "projecting..." << "\n";
      }
      Eigen::VectorXd sqrD;
      Eigen::MatrixX3d V_projected;
      igl::point_mesh_squared_distance(V_etc.leftCols(3).eval(), V0, F0, sqrD,
                                       Fi_containing_V_proj_out,
                                       V_projected); // Project
      V_etc.leftCols<3>() = V_projected;
    }

    if (verbose) {
      std::cout << "finished iter, faces amount: " << F_out.rows() << "\n";
    }
  }
  // return in Vattrs_out, removing the appended Vselection in V_etc
  Vattrs_out = V_etc(Eigen::all, Eigen::seq(0, Eigen::last - Eigen::fix<1>));
  // return in Vselection_out
  Vselection_out = (V_etc.rightCols<1>().array() >= selection_threshold)
                       .cast<int>()
                       .matrix();
}

// helper for finding barycentric hits with Fi_containing_V_proj
void barycentric_interp_on_Fi_containing_V_proj(
    const Eigen::MatrixXd &Vattrs_orig_in, const Eigen::MatrixX3i &F_orig_in,
    const Eigen::MatrixX3d &V_proj_in,
    const Eigen::VectorXi &Fi_containing_V_proj_in,
    Eigen::MatrixX3d &barycentric_coordinates_out,
    Eigen::MatrixXd &Vattrs_barylerped_out) {
  igl::barycentric_coordinates(
      V_proj_in,
      Vattrs_orig_in.leftCols<3>()(F_orig_in(Fi_containing_V_proj_in, 0),
                                   Eigen::all),
      Vattrs_orig_in.leftCols<3>()(F_orig_in(Fi_containing_V_proj_in, 1),
                                   Eigen::all),
      Vattrs_orig_in.leftCols<3>()(F_orig_in(Fi_containing_V_proj_in, 2),
                                   Eigen::all),
      barycentric_coordinates_out);
  igl::barycentric_interpolation(
      Vattrs_orig_in, F_orig_in, barycentric_coordinates_out,
      Fi_containing_V_proj_in, Vattrs_barylerped_out);
}