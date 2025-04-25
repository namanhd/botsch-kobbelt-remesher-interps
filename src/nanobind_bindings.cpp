#include "adaptive_remesh_botsch.h"
#include <Eigen/Core>
#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/tuple.h>

namespace nb = nanobind;

std::tuple<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::VectorXi>
remesh_botsch_with_interps(const nb::DRef<const Eigen::MatrixXd> &Vattrs,
                           const nb::DRef<const Eigen::MatrixXi> &F,
                           const nb::DRef<const Eigen::VectorXi> &Vselection,
                           double targetlen, double selection_threshold,
                           int iterations, bool project, int verbose) {
  Eigen::MatrixXd Vattrs_remesh;
  Eigen::MatrixXi F_remesh;
  Eigen::VectorXi Vselection_remesh;
  remesh_botsch(Vattrs, F, Vselection, targetlen, selection_threshold,
                iterations, project, verbose, Vattrs_remesh, F_remesh,
                Vselection_remesh);
  return std::make_tuple(Vattrs_remesh, F_remesh, Vselection_remesh);
}

std::tuple<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::VectorXi>
remesh_botsch_with_interps(const nb::DRef<const Eigen::MatrixXd> &Vattrs,
                           const nb::DRef<const Eigen::MatrixXi> &F,
                           const nb::DRef<const Eigen::VectorXi> &Vselection,
                           const nb::DRef<const Eigen::VectorXd> &targetlen,
                           double selection_threshold, int iterations,
                           bool project, int verbose) {
  Eigen::MatrixXd Vattrs_remesh;
  Eigen::MatrixXi F_remesh;
  Eigen::VectorXi Vselection_remesh;
  remesh_botsch(Vattrs, F, Vselection, targetlen, selection_threshold,
                iterations, project, verbose, Vattrs_remesh, F_remesh,
                Vselection_remesh);
  return std::make_tuple(Vattrs_remesh, F_remesh, Vselection_remesh);
}

std::tuple<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::VectorXi>
remesh_botsch_adaptive_with_interps(
    const nb::DRef<const Eigen::MatrixXd> &Vattrs,
    const nb::DRef<const Eigen::MatrixXi> &F,
    const nb::DRef<const Eigen::VectorXi> &Vselection, double epsilon,
    bool adaptive, double selection_threshold, int iterations, bool project,
    int verbose) {
  Eigen::MatrixXd Vattrs_remesh;
  Eigen::MatrixXi F_remesh;
  Eigen::VectorXi Vselection_remesh;
  remesh_botsch_adaptive(Vattrs, F, Vselection, epsilon, adaptive,
                         selection_threshold, iterations, project, verbose,
                         Vattrs_remesh, F_remesh, Vselection_remesh);
  return std::make_tuple(Vattrs_remesh, F_remesh, Vselection_remesh);
}

NB_MODULE(bkremeshlerps, m) {
  m.def("remesh_botsch_with_interps",
        (std::tuple<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::VectorXi>(*)(
            const nb::DRef<const Eigen::MatrixXd> &,
            const nb::DRef<const Eigen::MatrixXi> &,
            const nb::DRef<const Eigen::VectorXi> &, double, double, int, bool,
            int)) &
            remesh_botsch_with_interps,
        nb::arg("Vattrs"), nb::arg("F"), nb::arg("Vselection"),
        nb::arg("targetlen"), nb::arg("selection_threshold"),
        nb::arg("iterations"), nb::arg("project"), nb::arg("verbose") = 0,
        "Isotropic remeshing, with per-vertex target edge length. "
        "Features ability to select vertices allowed for "
        "remeshing and interpolate per-vertex attributes through remesh "
        "operations. "
        "(Overload allowing a constant target edge length for the whole mesh)");

  m.def("remesh_botsch_with_interps",
        (std::tuple<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::VectorXi>(*)(
            const nb::DRef<const Eigen::MatrixXd> &,
            const nb::DRef<const Eigen::MatrixXi> &,
            const nb::DRef<const Eigen::VectorXi> &,
            const nb::DRef<const Eigen::VectorXd> &, double, int, bool, int)) &
            remesh_botsch_with_interps,
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
        "computed via "
        "curvature.");
}
