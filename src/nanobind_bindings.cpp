#include "adaptive_remesh_botsch.h"
#include <Eigen/Core>
#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/tuple.h>

namespace nb = nanobind;

std::tuple<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::VectorXi, Eigen::VectorXi>
remesh_botsch_with_interps__consttargetlen(
    const nb::DRef<const Eigen::MatrixXd> &Vattrs,
    const nb::DRef<const Eigen::MatrixXi> &F,
    const nb::DRef<const Eigen::VectorXi> &Vselection, double targetlen,
    double selection_threshold, int iterations, bool project, int verbose) {
  Eigen::MatrixXd Vattrs_remesh;
  Eigen::MatrixX3i F_remesh;
  Eigen::VectorXi Vselection_remesh;
  Eigen::VectorXi Fi_containing_V_proj;
  remesh_botsch(Vattrs, F, Vselection, targetlen, selection_threshold,
                iterations, project, verbose, Vattrs_remesh, F_remesh,
                Vselection_remesh, Fi_containing_V_proj);
  return std::make_tuple(Vattrs_remesh, F_remesh, Vselection_remesh,
                         Fi_containing_V_proj);
}

std::tuple<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::VectorXi, Eigen::VectorXi>
remesh_botsch_with_interps__arrtarglen(
    const nb::DRef<const Eigen::MatrixXd> &Vattrs,
    const nb::DRef<const Eigen::MatrixXi> &F,
    const nb::DRef<const Eigen::VectorXi> &Vselection,
    const nb::DRef<const Eigen::VectorXd> &targetlen,
    double selection_threshold, int iterations, bool project, int verbose) {
  Eigen::MatrixXd Vattrs_remesh;
  Eigen::MatrixX3i F_remesh;
  Eigen::VectorXi Vselection_remesh;
  Eigen::VectorXi Fi_containing_V_proj;
  remesh_botsch(Vattrs, F, Vselection, targetlen, selection_threshold,
                iterations, project, verbose, Vattrs_remesh, F_remesh,
                Vselection_remesh, Fi_containing_V_proj);
  return std::make_tuple(Vattrs_remesh, F_remesh, Vselection_remesh,
                         Fi_containing_V_proj);
}

std::tuple<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::VectorXi, Eigen::VectorXi>
remesh_botsch_adaptive_with_interps(
    const nb::DRef<const Eigen::MatrixXd> &Vattrs,
    const nb::DRef<const Eigen::MatrixXi> &F,
    const nb::DRef<const Eigen::VectorXi> &Vselection, double epsilon,
    bool adaptive, double selection_threshold, int iterations, bool project,
    int verbose) {
  Eigen::MatrixXd Vattrs_remesh;
  Eigen::MatrixX3i F_remesh;
  Eigen::VectorXi Vselection_remesh;
  Eigen::VectorXi Fi_containing_V_proj;
  remesh_botsch_adaptive(Vattrs, F, Vselection, epsilon, adaptive,
                         selection_threshold, iterations, project, verbose,
                         Vattrs_remesh, F_remesh, Vselection_remesh,
                         Fi_containing_V_proj);
  return std::make_tuple(Vattrs_remesh, F_remesh, Vselection_remesh,
                         Fi_containing_V_proj);
}

std::tuple<Eigen::MatrixX3d, Eigen::MatrixXd>
barycentric_interp_on_Fi_containing_V_proj__binding(
    const nb::DRef<const Eigen::MatrixXd> &Vattrs_orig,
    const nb::DRef<const Eigen::MatrixX3i> &F_orig,
    const nb::DRef<const Eigen::MatrixX3d> &V_proj,
    const nb::DRef<const Eigen::VectorXi> &Fi_containing_V_proj) {
  Eigen::MatrixX3d barycentric_coordinates;
  Eigen::MatrixXd Vattrs_barylerped;
  barycentric_interp_on_Fi_containing_V_proj(
      Vattrs_orig, F_orig, V_proj, Fi_containing_V_proj,
      barycentric_coordinates, Vattrs_barylerped);
  return std::make_tuple(barycentric_coordinates, Vattrs_barylerped);
}

NB_MODULE(bkremeshlerps, m) {
  m.def("remesh_botsch_with_interps",
        &remesh_botsch_with_interps__consttargetlen, nb::arg("Vattrs"),
        nb::arg("F"), nb::arg("Vselection"), nb::arg("targetlen"),
        nb::arg("selection_threshold"), nb::arg("iterations"),
        nb::arg("project"), nb::arg("verbose") = 0,
        "Isotropic remeshing, with per-vertex target edge length. "
        "Features ability to select vertices allowed for "
        "remeshing and interpolate per-vertex attributes through remesh "
        "operations. "
        "(Overload allowing a constant target edge length for the whole mesh)");

  m.def("remesh_botsch_with_interps", &remesh_botsch_with_interps__arrtarglen,
        nb::arg("Vattrs"), nb::arg("F"), nb::arg("Vselection"),
        nb::arg("targetlen"), nb::arg("selection_threshold"),
        nb::arg("iterations"), nb::arg("project"), nb::arg("verbose") = 0,
        "Isotropic remeshing, with per-vertex target edge length. "
        "Features ability to select vertices allowed for "
        "remeshing and interpolate per-vertex attributes through remesh "
        "operations.");

  m.def("remesh_botsch_adaptive_with_interps",
        &remesh_botsch_adaptive_with_interps, nb::arg("Vattrs"), nb::arg("F"),
        nb::arg("Vselection"), nb::arg("epsilon"), nb::arg("adaptive"),
        nb::arg("selection_threshold"), nb::arg("iterations"),
        nb::arg("project"), nb::arg("verbose") = 0,
        "Adaptive isotropic remeshing, with target edge length/sizing field "
        "computed via curvature.");

  m.def(
      "barycentric_interp_on_Fi_containing_V_proj",
      &barycentric_interp_on_Fi_containing_V_proj__binding,
      nb::arg("Vattrs_orig"), nb::arg("F_orig"), nb::arg("V_proj"),
      nb::arg("Fi_containing_V_proj"),
      "Barycentric interpolate original per-vertex attributes onto the "
      "remeshed projected vertices V_proj given the faces in F_orig that "
      "contain V_proj. This is another way to obtain interpolated attributes "
      "through the remesh (usable after the remesh_* functions), though with "
      "less guaranteed numerical stability since we use barycentric coords of "
      "the new remeshed vertices' closest faces on the source mesh rather than "
      "interpolation that happens continuously at each remeshing operation");
}
