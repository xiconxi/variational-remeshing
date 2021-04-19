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
    Eigen::SparseMatrix<double> mass1, mass2;
    igl::massmatrix(V_, F_, igl::MASSMATRIX_TYPE_VORONOI, mass1);
    igl::massmatrix(UV_, F_, igl::MASSMATRIX_TYPE_VORONOI, mass2);
    density_ = mass1.diagonal().array()/mass2.diagonal().array();
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

std::pair<Eigen::Vector2d, double> density_integral(Eigen::MatrixX2d& poly, Eigen::VectorXd& d) {
    double mean_0 = poly.col(0).dot(d)/d.mean();
    double mean_1 = poly.col(1).dot(d)/d.mean();
    double area = 0;
    Eigen::Vector3d v0 , v1;
    v0.setZero();
    v1.setZero();
    for(int i = 1; i < poly.size() -1; i++) {
        v0.topRows(2) = poly.row(i) - poly.row(0);
        v1.topRows(2) = poly.row(i+1) - poly.row(0);
        area += v0.cross(v1).norm();
    }

    return {Eigen::Vector2d(mean_0, mean_1), area};
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
    for(int i = 0; i < diagram.numsites; i++) {
        clip.clear();
        for(auto e = sites[i].edges; e; e = e->next)
            clip.push_back(e->pos[0].x+e->pos[0].y*1.0i);

        Eigen::Vector2d centroid, centroid_sum;
        double weight = 0, weight_sum = 0;

        Eigen::VectorXd sqrD;
        Eigen::VectorXi FI;
        Eigen::MatrixXd C, P;

        centroid_sum.setZero();
        for(int fi = 0; fi < F_.rows(); fi++) {
            subj[0] = V[F_(fi, 0)];
            subj[1] = V[F_(fi, 1)];
            subj[2] = V[F_(fi, 2)];
            // clip the triangle with voronoi cell: subj&clip -> solution
            clipping_sutherland_hodgman(subj, clip, solution);
            if(solution.size() == 0) continue;

            P.resize(solution.size(), 2);
            for(int s_i = 0; s_i < P.rows(); s_i++) {
                P(s_i, 0) = solution[s_i].real();
                P(s_i, 1) = solution[s_i].imag();
            }

            //compute the barycentric coord of the voronoi sites -> C
            tree.squared_distance(UV_, F_, P, sqrD, FI, C);
            Eigen::VectorXd d(FI.size());
            for(int k = 0; k < d.size(); k++)
                d[k] = density(FI[k], C.row(k));




            std::cout << sqrD.transpose() << std::endl;
        }
    }


    jcv_diagram_free( &diagram );

}