//
// Created by pupa on 4/15/21.
//

#include "weighted_cvd.h"

#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/lscm.h>
#include <igl/AABB.h>

#include <igl/massmatrix.h>
#include "polygon_clip.h"

#define JCV_REAL_TYPE double
#define JC_VORONOI_IMPLEMENTATION
#include "jc_voronoi.h"



void VariationalRemesher::lscm_parameterization() {
    // Find the open boundary
    Eigen::VectorXi bnd;
    igl::boundary_loop(F_, bnd);

    // Map the boundary to a circle, preserving edge proportions
    Eigen::MatrixXd bnd_uv;
    igl::map_vertices_to_circle(V_, bnd,bnd_uv);

    igl::lscm(V_, F_, bnd, bnd_uv, UV_);
}


void VariationalRemesher::compute_local_density() {
    Eigen::SparseMatrix<double> mass;
    igl::massmatrix(V_, F_, igl::MASSMATRIX_TYPE_VORONOI, mass);
    density_ = mass.diagonal();
}

void VariationalRemesher::initial_sampling(size_t n_sample) {
    samples_.resize(n_sample, 2);
    for (int i = 0; i < n_sample; i++) {
        double r = std::sqrt((rand() % 1000) / 1000.0);
        double angle = (rand() % 1000) / 1000.0 * M_PI * 2;
        samples_(i, 0) = r * std::cos(angle);
        samples_(i, 1) = r * std::sin(angle);
    }
}

double VariationalRemesher::density_integral(Eigen::MatrixX2d& poly, Eigen::VectorXd d) {

}

void VariationalRemesher::lloyd_relaxation(size_t n_iters) {
    Eigen::VectorXcd V(V_.rows());
    for(int i = 0; i < V_.rows(); i++)
        V[i] = V_(i, 0) + V_(i, 1)*1i;

    std::vector<jcv_point> points(samples_.rows());
    for(int i = 0; i < samples_.rows(); i++) {
        points[i].x = samples_(i, 0);
        points[i].y = samples_(i, 1);
    }
    jcv_diagram diagram;
    memset(&diagram, 0, sizeof(jcv_diagram));
    jcv_diagram_generate(points.size(), points.data(), 0, 0, &diagram );

    const jcv_site* sites = jcv_diagram_get_sites(&diagram);
    std::vector<Complex> subj(3), clip, solution;
    clip.reserve(10);
    igl::AABB<Eigen::MatrixXd, 2> tree;
    tree.init(UV_, F_);
    Eigen::VectorXd sqrD;
    Eigen::VectorXi I;
    Eigen::MatrixXd C, P;
    for(int i = 0; i < diagram.numsites; i++) {
        clip.clear();
        for(auto e = sites[i].edges; e; e = e->next)
            clip.push_back(e->pos[0].x+e->pos[0].y*1.0i);
        for(int fi = 0; fi < F_.rows(); fi++) {
            subj[0] = V[F_(fi, 0)];
            subj[1] = V[F_(fi, 1)];
            subj[2] = V[F_(fi, 2)];
            clipping_sutherland_hodgman(subj, clip, solution);
            if(solution.size() == 0) continue;
            P.resize(solution.size(), 2);
            for(int s_i = 0; s_i < P.rows(); s_i++) {
                P(s_i, 0) = solution[s_i].real();
                P(s_i, 1) = solution[s_i].imag();
            }
            tree.squared_distance(UV_, F_, P, sqrD, I, C);
            std::cout << sqrD.transpose() << std::endl;
        }
    }


    jcv_diagram_free( &diagram );

}