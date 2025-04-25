#include <igl/per_vertex_normals.h>
// #include <igl/principal_curvature.h>
// #include <igl/avg_edge_length.h>
// #include <igl/massmatrix.h>
#include <igl/adjacency_list.h>
#include <igl/barycenter.h>
// #include <igl/pinv.h>
// #include <igl/writeOBJ.h>
// #include <igl/edges.h>
// #include <Eigen/SparseCore>
// #include <igl/adjacency_matrix.h>
#include <igl/per_face_normals.h>
// #include <igl/unique_edge_map.h>
#include <igl/vertex_triangle_adjacency.h>
// #include <igl/flip_edge.h>
#include <igl/remove_duplicate_vertices.h>


void tangential_smoothing(Eigen::MatrixXd &V_etc, Eigen::MatrixXi &F, double selthresh, Eigen::VectorXi &feature, Eigen::VectorXd &lambda, 
    const Eigen::VectorXd &sizingField, bool use_sizing_field_weighted_avg)
{
    using namespace Eigen;
    Eigen::MatrixXd N;
    Eigen::VectorXd dblA;
    std::vector<std::vector<int>> A;
    Matrix3d I, NN;
    I.setIdentity();
    Eigen::MatrixXd SV;
    Eigen::MatrixXi SVI, SVJ;

    int n = V_etc.rows();  
    const int n_attribs = V_etc.cols();  
    int m = F.rows();

    // igl::doublearea(V,F,dblA);

    // std::vector<double> vertex_areas;
    // vertex_areas.setZero(m);

    // for (int j = 0; j < m; j++) {
    //     vertex_areas[F(j,0)] = vertex_areas[F(j,0)] + (abs(dblA(j))/6);
    //     vertex_areas[F(j,1)] = vertex_areas[F(j,1)] + (abs(dblA(j))/6);
    //     vertex_areas[F(j,2)] = vertex_areas[F(j,2)] + (abs(dblA(j))/6);
    // }

    Eigen::MatrixXd N_before, N_after;
    igl::adjacency_list(F, A);

    int num_feat = feature.size();
    std::vector<bool> is_feature_vertex;
    is_feature_vertex.resize(n);

    for (int s = 0; s < num_feat; s++)
    {
        is_feature_vertex[feature(s)] = true;
    }

    // Q.resize(n, 3);
    // P.resize(n, 3);
    //           Eigen::MatrixXd N;
    igl::per_vertex_normals(V_etc.leftCols<3>(), F, N);

    std::vector<std::vector<int>> V2F, _1;
    igl::vertex_triangle_adjacency(n, F, V2F, _1);
    Eigen::VectorXd area;
    Eigen::MatrixXd BC;
    igl::doublearea(V_etc.leftCols<3>(), F, area);
    area /= 2;
    igl::barycenter(V_etc, F, BC);  // n-d barycenter, not just 3D!

    Eigen::RowVectorXd q(n_attribs), p(n_attribs);
    MatrixXd NN_(n_attribs, n_attribs);
    NN_.setIdentity();
    for (int i = 0; i < n; i++)
    {
        // skip if vertex is not selected
        // std::cout << "i: " << i << "/" << n << "\n";
        // bool is_feature = is_feature_vertex[i];
        bool allow_this_vertex = V_etc(i, Eigen::last) >= selthresh && (!is_feature_vertex[i]);
        if (allow_this_vertex)
        {
            double denom = 0;
            q.setZero();
            p.setZero();
            if (use_sizing_field_weighted_avg) {
                /* yotam's */
                RowVectorXd nom(n_attribs);
                nom.setZero();
                for (int j = 0; j < V2F[i].size(); j++)
                {
                    double cur_denom = (area(V2F[i][j]) * (sizingField(F(V2F[i][j], 0)) + sizingField(F(V2F[i][j], 1)) + sizingField(F(V2F[i][j], 2))) / 3);
                    denom = denom + cur_denom;
                    nom = nom + (cur_denom * BC.row(V2F[i][j]));
                }

                // so this is 
                // sum_{f in adjacentFaces}((area of f * mean of sizing of the verts of face f) * barycenter of face f)
                // divided by (sum_{f in adjacentFaces}(area of f * mean of sizing of the verts of face f))
                //  i.e. q is a weighted average of the barycenters of adjacent faces
                // (where the weight of a face is its area times its mean sizing field value)
                
                // std::cout << "V[0]: " << V.row(0) << "\n";
                // std::cout << "nom: " << nom << "\n";
                // std::cout << "denom: " << denom << "\n";
                q = nom.array() / denom;
                q = q.matrix();
            }
            else {
                /* silvia's */
                for (int j = 0; j < A[i].size(); j++) {
                    q = q + (V_etc.row(A[i][j]) / A[i].size());
                    // q = q + (V.row(A[i][j])*vertex_areas[A[i][j]]);
                    // std::cout << q << std::endl;
                    // denominator = denominator + vertex_areas[A[i][j]];
                }
            }

            // namanh NOTE: N is (nverts, 3), N.row(i) is RowVector3d, transpose gets col Vector3d, then mult with row Vector 3d
            //  (N.row(i).transpose() * (N.row(i)))  is a (1,3).T @ (1,3) which is (3,1)@(1,3) which is (3,3)
            // so this is the outer product of the normal with itself...
            
            NN = lambda(i) * (Eigen::MatrixXd::Identity(3, 3) - N.row(i).transpose() * (N.row(i)));
            NN_.topLeftCorner<3, 3>() = NN;

            
            // then you take this NN (3,3) @ (verts[i] - q[i])
            // this result is always 
            // (v_etc[i]-q) minus ((vnormals[i] dot (v_etc[i]-q)) .* (v_etc[i]-q))
            // and then 
            // p = v_etc[i] - (v_etc[i]-q - ((vnormals[i] dot (v_etc[i]-q)) .* (v_etc[i]-q)))
            // p = q + ((vnormals[i] dot (v_etc[i]-q)) .* (v_etc[i]-q))
            // lol what..

            p = (V_etc.row(i).transpose() - (NN_ * (V_etc.row(i).transpose() - q.transpose()))).transpose();
            V_etc.row(i) = p;

            // p = q;
            // std::cout << N.row(i) << std::endl;
            // Eigen::RowVectorXd padded_normal(n_attribs);
            // padded_normal.setConstant(1.0);
            // padded_normal.head<3>() = N.row(i);
            // // double ndotq = N.row(i) * q;
            // V_etc.row(i) = q + (padded_normal * q) * (V_etc.row(i) - q);
            // std::cout << "q: " << q << "\n";
            // std::cout << "V_etc.row(i): " << V_etc.row(i) << "\n";
            // std::cout << p << std::endl;

            // igl::per_face_normals(V_projected,F,Eigen::Vector3d(0,0,0),N_after);
            //            for (int j = 0; j < m ; j++) {
            //                if (N_before.row(j).dot(N_after.row(j)) < 0) {
            //                    // std::cout << "Avoided face flipping, I think." << std::endl;
            //                    V.row(i) = V_fixed.row(i);
            //                }
            //            }
        }
    }
    //        igl::remove_duplicate_vertices(V,0,SV,SVI,SVJ);
    //        std::cout << V.rows()-SV.rows() << std::endl;
    //	igl::writeOBJ("pre-project.obj",V,F);
    //
    
    //	igl::writeOBJ("post-project.obj",V,F);
    //    igl::remove_duplicate_vertices(V,0,SV,SVI,SVJ);
    //    std::cout << V.rows()-SV.rows() << std::endl;
    //   std::cout << "not projecting!" << std::endl;
}