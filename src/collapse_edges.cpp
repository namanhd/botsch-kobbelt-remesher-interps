// #include <igl/per_vertex_normals.h>
// #include <igl/principal_curvature.h>
// #include <igl/avg_edge_length.h>
// #include <igl/massmatrix.h>
#include <Eigen/Core>
#include <igl/adjacency_list.h>
// #include <igl/per_face_normals.h>
// #include <igl/barycenter.h>
// #include <igl/pinv.h>
// #include <igl/edges.h>
// #include <Eigen/SparseCore>
// #include <igl/adjacency_matrix.h>
// #include <igl/edge_flaps.h>
// #include <igl/unique_edge_map.h>
// #include <igl/vertex_triangle_adjacency.h>
// #include <igl/collapse_edge.h>
#include <igl/is_edge_manifold.h>
// #include <igl/edge_collapse_is_valid.h>
// #include <igl/C_STR.h>
// #include <igl/circulation.h>
#include <cassert>
#include <igl/decimate.h>
#include <igl/decimate_callback_types.h>
#include <igl/infinite_cost_stopping_condition.h>
#include <igl/remove_unreferenced.h>
#include <igl/shortest_edge_and_midpoint.h>

template <typename DerivedV_etc, typename DerivedF, typename DerivedFeature,
          typename DerivedHigh, typename DerivedLow>
void collapse_edges(Eigen::MatrixBase<DerivedV_etc> &V_etc,
                    Eigen::MatrixBase<DerivedF> &F, double selthresh,
                    Eigen::MatrixBase<DerivedFeature> &feature,
                    Eigen::MatrixBase<DerivedHigh> &high,
                    Eigen::MatrixBase<DerivedLow> &low) {
  Eigen::MatrixXi E, uE, EI, EF;
  Eigen::VectorXi EMAP, I, J;
  // Eigen::VectorXd data;
  Eigen::Matrix<typename DerivedV_etc::Scalar, Eigen::Dynamic, Eigen::Dynamic>
      U;
  Eigen::MatrixXi G;
  // namanh: note that V_etc is shape (n_verts, 3+1+attribs), contains
  // coordinates as well as other features being interpolated (1 for the
  // selection "bool", attribs for the rest) so we should recover V with only
  // 3D coords (first 3 coords) for geometric ops
  // seems like in Eigen the assignment operator = does a copy so it's better to
  // just slice every time we need a [:, :3] i.e.  (Eigen::all, Eigen::seq(0,3))
  // or leftCols<3>() which is actually zero-overhead
  std::vector<std::vector<int>> uE2E;
  std::vector<std::vector<int>> vertex_face_adjacency;
  std::vector<int> small_edges;
  int e1, e2, f1, f2, e;
  int n = V_etc.rows();
  const int n_attribs_minus1 = V_etc.cols() - 1;

  int num_feature = feature.size();
  std::vector<std::vector<int>> A;
  igl::adjacency_list(F, A);

  std::vector<bool> is_feature_vertex;
  is_feature_vertex.resize(n);

  for (int s = 0; s < num_feature; s++) {
    is_feature_vertex[feature(s)] = true;
  }

  // igl::is_edge_manifold(F);

  igl::decimate_stopping_condition_callback stopping_condition;

  // these are the original vertex indices of the edge to collapse
  // to be filled by pre_collapse_lambda to be then used by post_collapse_lambda
  // since the E array that post_collapse_lambda has access to has been mangled
  // by igl's collapse; we must store our own "original edge's vertex indices"
  // to fix up our own V_etc
  int e_v1 = -1, e_v2 = -1;
  // criterion for collapse can be || and not && selected, but we should
  // probably only set the result to be selected if both were selected
  bool both_selected = false;

  igl::decimate_pre_collapse_callback pre_collapse_lambda =
      [&V_etc, selthresh, &e_v1, &e_v2, &both_selected](
          const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
          const Eigen::MatrixXi &E, const Eigen::VectorXi &EMAP,
          const Eigen::MatrixXi &EF, const Eigen::MatrixXi &EI,
          const igl::min_heap<std::tuple<double, int, int>> &Q,
          const Eigen::VectorXi &EQ, const Eigen::MatrixXd &C, const int e) {
        bool sel1 = (V_etc(E(e, 0), Eigen::last) >= selthresh);
        bool sel2 = (V_etc(E(e, 1), Eigen::last) >= selthresh);
        bool selected = sel1 && sel2;
        bool both_selected = sel1 && sel2;
        if (selected) {
          e_v1 = E(e, 0);
          e_v2 = E(e, 1);
        }
        return selected;
      };

  igl::decimate_post_collapse_callback post_collapse_lambda =
      [&V_etc, selthresh, n_attribs_minus1, &e_v1, &e_v2, &high, &low,
       &both_selected](const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                       const Eigen::MatrixXi &E, const Eigen::VectorXi &EMAP,
                       const Eigen::MatrixXi &EF, const Eigen::MatrixXi &EI,
                       const igl::min_heap<std::tuple<double, int, int>> &Q,
                       const Eigen::VectorXi &EQ, const Eigen::MatrixXd &C,
                       const int e, const int e1, const int e2, const int f1,
                       const int f2, const bool collapsed) {
        if (collapsed) {
          assert(e_v1 >= 0 && e_v2 >= 0);
          // do what igl collapse_edge.cpp does but on V_etc (igl::decimate
          // and igl::collapse_edge only operated on V) we need to do the same
          // things to V_etc. But fortunately that's not a lot of things
          Eigen::RowVectorXd p_etc = (V_etc.row(e_v1) + V_etc.row(e_v2)) / 2;
          // if (both_selected)
          // {
          //     p_etc(n_attribs_minus1) = 1;
          // }
          V_etc.row(e_v1) = p_etc;
          V_etc.row(e_v2) = p_etc;
          // that's it, then afterwards we can use the I array returned from
          // igl's decimate call (which is the J from remove_unreferenced) on V;
          // this I has the indices of the surviving vertices to extract out
          // of both V and V_etc
          // wait we should do this with high and low too
          double p_high = (high(e_v1) + high(e_v2)) / 2;
          double p_low = (low(e_v1) + low(e_v2)) / 2;
          high(e_v1) = p_high;
          high(e_v2) = p_high;
          low(e_v1) = p_low;
          low(e_v2) = p_low;
          e_v1 = -1;
          e_v2 = -1;
        }
      };

  igl::decimate_cost_and_placement_callback shortest_edge_and_midpoint_lambda =
      [&A, &feature, &low, &high, &is_feature_vertex](
          const int e, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
          const Eigen::MatrixXi &E, const Eigen::VectorXi &EMAP,
          const Eigen::MatrixXi &EF, const Eigen::MatrixXi &EI, double &cost,
          Eigen::RowVectorXd &p) {
        igl::shortest_edge_and_midpoint(e, V, F, E, EMAP, EF, EI, cost, p);
        if (is_feature_vertex[E(e, 0)] || is_feature_vertex[E(e, 1)]) {
          cost = std::numeric_limits<double>::infinity();
          return;
        }
        if ((V.row(E(e, 0)) - V.row(E(e, 1))).norm() >
            ((low(E(e, 0)) + low(E(e, 1))) / 2)) {
          cost = std::numeric_limits<double>::infinity();
          return;
        }
        for (int i = 0; i < A[E(e, 1)].size(); i++) {
          if ((V.row(A[E(e, 1)][i]) - p).norm() > high(E(e, 1))) {
            cost = std::numeric_limits<double>::infinity();
            return;
          }
        }
        for (int r = 0; r < A[E(e, 0)].size(); r++) {
          if ((V.row(A[E(e, 0)][r]) - p).norm() > high(E(e, 0))) {
            cost = std::numeric_limits<double>::infinity();
            return;
          }
        }
        // consider each face

        // BUILD N
        int ccw = E(e, 0) > E(e, 1);
        std::vector<int> N;
        N.reserve(6);
        const int m = EMAP.size() / 3;
        assert(m * 3 == EMAP.size());
        const auto &step = [&](const int e, const int ff, int &ne, int &nf) {
          assert((EF(e, 1) == ff || EF(e, 0) == ff) && "e should touch ff");
          // const int fside = EF(e,1)==ff?1:0;
          const int nside = EF(e, 0) == ff ? 1 : 0;
          const int nv = EI(e, nside);
          // get next face
          nf = EF(e, nside);
          // get next edge
          const int dir = ccw ? -1 : 1;
          ne = EMAP(nf + m * ((nv + dir + 3) % 3));
        };
        // Always start with first face (ccw in step will be sure to turn right
        // direction)
        const int f0 = EF(e, 0);
        int fi = f0;
        int ei = e;
        while (true) {
          step(ei, fi, ei, fi);
          N.push_back(fi);
          // back to start?
          if (fi == f0) {
            break;
          }
        }
        for (const int f : N) {
          // std::cout << f << std::endl;
          if (f == 0 || f == N.size() - 1) {

            // skip
            continue;
          }
          // Grab the three corners of the face
          Eigen::RowVector3d p_before[3], p_after[3];
          for (int c = 0; c < 3; c++) {
            // vertex index
            //                    std::cout << e << std::endl;
            //                    std::cout << f << std::endl;
            //                    std::cout << c << std::endl;
            const int v = F(f, c);
            if (v == E(e, 0) || v == E(e, 1)) {
              p_after[c] = p;
            } else {
              p_after[c] = V.row(v);
            }
            p_before[c] = V.row(v);
          }
          const Eigen::RowVector3d n_before =
              ((p_before[1] - p_before[0]).cross(p_before[2] - p_before[0]))
                  .normalized();
          const Eigen::RowVector3d n_after =
              ((p_after[1] - p_after[0]).cross(p_after[2] - p_after[0]))
                  .normalized();
          if (n_before.dot(n_after) < n_after.norm() / 2) {
            cost = std::numeric_limits<double>::infinity();
          }
        }

        // std::cout << "Mathed!" << std::endl;
      };

  igl::infinite_cost_stopping_condition(shortest_edge_and_midpoint_lambda,
                                        stopping_condition);

  // std::cout << "??" << std::endl;
  igl::decimate(V_etc.template leftCols<3>().eval(), F,
                shortest_edge_and_midpoint_lambda, stopping_condition,
                pre_collapse_lambda, post_collapse_lambda, U, G, J, I);
  // std::cout << "!!" << std::endl;

  // Eigen::VectorXd high_new, low_new;
  Eigen::VectorXi feature_new;
  feature_new.resize(num_feature);
  // high_new.resize(U.rows());
  // low_new.resize(U.rows());
  int j = 0;
  for (int s = 0; s < U.rows(); s++) {
    // high_new(s) = high(I(s));
    // low_new(s) = low(I(s));
    if (is_feature_vertex[I(s)]) {
      feature_new(j) = s;
      j = j + 1;
    }
  }

  // PLACEHOLDER

  V_etc = V_etc(I, Eigen::all).eval();
  F = G;
  high = high(I, Eigen::all).eval();
  low = low(I, Eigen::all).eval();
  // high = high_new;
  // low = low_new;
  feature = feature_new;
}

// g++ -I/usr/local/libigl/external/eigen -I/usr/local/libigl/include -std=c++11
// -framework Accelerate main.cpp remesh_botsch.cpp -o main
