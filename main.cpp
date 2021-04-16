#include <igl/boundary_loop.h>
#include <igl/harmonic.h>

#include <igl/knn.h>
#include <igl/intersect.h>
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/lscm.h>
#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>

#include "weighted_cvd.h"



Eigen::MatrixXd V;
Eigen::MatrixX3i F;
Eigen::MatrixXd V_uv;

bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier) {
    if (key == '1'){
        // Plot the 3D mesh
        viewer.data().set_mesh(V,F);
        viewer.core().align_camera_center(V,F);
    }
    else if (key == '2') {
        // Plot the mesh in 2D using the UV coordinates as vertex coordinates
        viewer.data().set_mesh(V_uv,F);
        viewer.core().align_camera_center(V_uv,F);
    }

    viewer.data().compute_normals();

    return false;
}

#include "polygon_clip.h"


int main(int argc, char *argv[])
{
    // Load a mesh in OFF format
    igl::read_triangle_mesh(argv[1] , V, F);

    VariationalRemesher remesher(V, F);
    remesher.lscm_parameterization();
    remesher.compute_local_density();
    remesher.initial_sampling(100);
    remesher.lloyd_relaxation(1);

    return 0;

    // Find the open boundary
    Eigen::VectorXi bnd;
    igl::boundary_loop(F,bnd);

    // Map the boundary to a circle, preserving edge proportions
    Eigen::MatrixXd bnd_uv;
    igl::map_vertices_to_circle(V,bnd,bnd_uv);

    igl::lscm(V, F, bnd, bnd_uv, V_uv);
    // Harmonic parametrization for the internal vertices
//    igl::harmonic(V,F,bnd,bnd_uv,1,V_uv);

    // Scale UV to make the texture more clear
    V_uv *= 5;

    // Plot the mesh
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V, F);
    viewer.data().set_uv(V_uv);
    viewer.callback_key_down = &key_down;

    // Disable wireframe
    viewer.data().show_lines = false;

    // Draw checkerboard texture
    viewer.data().show_texture = true;

    // Launch the viewer
    viewer.launch();
}
