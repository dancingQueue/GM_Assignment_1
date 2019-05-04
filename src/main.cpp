#include <iostream>
#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/per_corner_normals.h>
/*** insert any libigl headers here ***/

using namespace std;
using Viewer = igl::opengl::glfw::Viewer;

// Vertex array, #V x3
Eigen::MatrixXd V;
// Face array, #F x3
Eigen::MatrixXi F;
// Per-face normal array, #F x3
Eigen::MatrixXd FN;
// Per-vertex normal array, #V x3
Eigen::MatrixXd VN;
// Per-corner normal array, (3#F) x3
Eigen::MatrixXd CN;
// Vectors of indices for adjacency relations
std::vector<std::vector<int> > VF, VFi, VV;
// Integer vector of component IDs per face, #F x1
Eigen::VectorXi cid;
// Per-face color array, #F x3
Eigen::MatrixXd component_colors_per_face;

bool callback_key_down(Viewer& viewer, unsigned char key, int modifiers) {
    if (key == '1') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        // Add your code for computing vertex to face relations here;
        // store in VF,VFi.
        igl::vertex_triangle_adjacency(V.rows(),F,VF,VFi);

        for (auto row : VF) {
            cout << "\n";
            for (auto elem : row) {
                cout << elem << " ";
            }
            cout << "\n";
        }

    }

    if (key == '2') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        // Add your code for computing vertex to vertex relations here:
        // store in VV.
        igl::adjacency_list(F, VV);

        for (auto row : VV) {
            cout << "\n";
            for (auto elem : row) {
                cout << elem << " ";
            }
            cout << "\n";
        }
    }

    if (key == '3') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        FN.setZero(F.rows(),3);
        // Add your code for computing per-face normals here: store in FN.
        igl::per_face_normals(V,F,FN);
        // Set the viewer normals.
        viewer.data().set_normals(FN);
    }

    if (key == '4') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        VN.setZero(V.rows(), 3);
        // Add your code for computing per-vertex normals here: store in VN.
        igl::per_vertex_normals(V,F,VN);
        // Set the viewer normals.
        viewer.data().set_normals(VN);
    }

    if (key == '5') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        // Add your code for computing per-corner normals here: store in CN.
        CN.setZero(3 * F.rows(), 3);
        igl::per_corner_normals(V, F, 20, CN);
        //Set the viewer normals
        viewer.data().set_normals(CN);
    }

    if (key == '6') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        component_colors_per_face.setZero(F.rows(),3);
        // Add your code for computing per-face connected components here:
        // store the component labels in cid.

        // Compute colors for the faces based on components, storing them in
        // component_colors_per_face.

        // Set the viewer colors
        viewer.data().set_colors(component_colors_per_face);
    }

    if (key == '7') {
        Eigen::MatrixXd Vout=V;
        Eigen::MatrixXi Fout=F;
        // Add your code for sqrt(3) subdivision here.

        // Set up the viewer to display the new mesh
        V = Vout; F = Fout;
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
    }

    // Loop subdivision scheme
    if (key == '8') {
        std::vector<std::vector<int> > VVt;
        igl::adjacency_list(F, VVt);

        VV = VVt;
        Eigen::MatrixXd Vout(4 * V.rows(), 3);
        Eigen::MatrixXi Fout(4 * F.rows(), 3);

        // Compute new "even" vertices
        for (int i = 0; i < V.rows(); i++) {
            if (i >= VV.size()) break;
            int valence = VV[i].size();

            double beta = 3. / (8 * valence);

            Vout.row(i) = (1 - beta * valence) * V.row(i);

            for (auto u : VV[i]) {
                Vout.row(i) += beta * V.row(u);
            }
        }
        // Compute new faces and "odd" vertices
        int lastIndex = V.rows();
        Eigen::SparseMatrix<int> edgeToV(lastIndex, lastIndex);

        for (int i = 0; i < F.rows(); i += 1) {
            auto A = V.row(F(i, 0));
            auto B = V.row(F(i, 1));
            auto C = V.row(F(i, 2));

            auto A_index = F(i, 0);
            auto B_index = F(i, 1);
            auto C_index = F(i, 2);

            auto midAB_index = lastIndex;

            if (edgeToV.coeff(A_index, B_index) > 0) {
                midAB_index = edgeToV.coeff(A_index, B_index);
            } else {
                 edgeToV.insert(A_index, B_index) = lastIndex;
                 edgeToV.insert(B_index, A_index) = lastIndex;
                 lastIndex += 1;
            }

            auto midBC_index = lastIndex;

            if (edgeToV.coeff(B_index, C_index) > 0) {
                midBC_index = edgeToV.coeff(B_index, C_index);
            } else {
                 edgeToV.insert(B_index, C_index) = lastIndex;
                 edgeToV.insert(C_index, B_index) = lastIndex;
                 lastIndex += 1;
            }

            auto midAC_index = lastIndex;

            if (edgeToV.coeff(A_index, C_index) > 0) {
                midAC_index = edgeToV.coeff(A_index, C_index);
            } else {
                 edgeToV.insert(A_index, C_index) = lastIndex;
                 edgeToV.insert(C_index, A_index) = lastIndex;
                 lastIndex += 1;
            }
            
            Vout.row(midAB_index) = 3. / 8. * (A + B);
            Vout.row(midBC_index) = 3. / 8. * (B + C);
            Vout.row(midAC_index) = 3. / 8. * (A + C);

            int count = 0;
            for (auto nU : VV[A_index]) {
                for (auto nV : VV[B_index]) {
                    if (count == 2) break;
                    if (nU == nV && (nU != A_index) && (nV != B_index)) {
                        Vout.row(midAB_index) += (1. / 8. * V.row(nU));
                        count += 1;
                    }
                }
                if (count == 2) break;
            }

            count = 0;
            for (auto nU : VV[B_index]) {
                for (auto nV : VV[C_index]) {
                    if (count == 2) break;
                    if (nU == nV && (nU != B_index) && (nV != C_index)) {
                        Vout.row(midBC_index) += (1. / 8. * V.row(nU));
                        count += 1;
                    }
                }
                if (count == 2) break;
            }

            count = 0;
            for (auto nU : VV[A_index]) {
                for (auto nV : VV[C_index]) {
                    if (count == 2) break;
                    if (nU == nV && (nU != A_index) && (nV != C_index)) {
                        Vout.row(midAC_index) += (1. / 8. * V.row(nU));
                        count += 1;
                    }

                }
                if (count == 2) break;
            }

            Fout(4 * i, 0) = midAB_index; 
            Fout(4 * i, 1) = midBC_index;
            Fout(4 * i, 2) = midAC_index;
            
            Fout(4 * i + 1, 1) = midAB_index; 
            Fout(4 * i + 1, 0) = F(i, 0);
            Fout(4 * i + 1, 2) = midAC_index;

            Fout(4 * i + 2, 0) = midAB_index; 
            Fout(4 * i + 2, 2) = midBC_index;
            Fout(4 * i + 2, 1) = F(i, 1);

            Fout(4 * i + 3, 2) = F(i, 2); 
            Fout(4 * i + 3, 1) = midBC_index;
            Fout(4 * i + 3, 0) = midAC_index;

        }



        // Set up the viewer to display the new mesh
        V = Vout; F = Fout;
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
    }

    return true;
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        cout << "Usage ex1_bin mesh.obj" << endl;
        //exit(0);
        V= (Eigen::MatrixXd(8,3)<<
    0.0,0.0,0.0,
    0.0,0.0,1.0,
    0.0,1.0,0.0,
    0.0,1.0,1.0,
    1.0,0.0,0.0,
    1.0,0.0,1.0,
    1.0,1.0,0.0,
    1.0,1.0,1.0).finished();
  F = (Eigen::MatrixXi(12,3)<<
    1,7,5,
    1,3,7,
    1,4,3,
    1,2,4,
    3,8,7,
    3,4,8,
    5,7,8,
    5,8,6,
    1,5,6,
    1,6,2,
    2,6,8,
    2,8,4).finished().array()-1;
    } else {
        // Read mesh
        igl::readOFF(argv[1],V,F);
    }


    // Show the mesh
    Viewer viewer;
    viewer.callback_key_down = callback_key_down;
    callback_key_down(viewer, '1', 0);

    viewer.launch();
}
