//
// Created by pupa on 4/15/21.
//

#pragma once

#include <Eigen/Dense>
#include <vector>

using std::vector;
using CCWPoly = vector<Eigen::RowVector2d>;

bool to_left_test(const Eigen::RowVector2d& p, const Eigen::RowVector2d& s, const Eigen::RowVector2d& t);

void convex_poly_intersect(const CCWPoly& subj, const CCWPoly& clip, CCWPoly& poly);

class VariationalRemesher{
public:
    explicit VariationalRemesher(Eigen::MatrixXd& V, Eigen::MatrixX3i& F)
            :V_(V), F_(F){}

    void lscm_parameterization();

    void compute_local_density();

    void initial_sampling(size_t n_sample);

    void lloyd_relaxation(size_t n_iters);

    double inline density(size_t fi, Eigen::RowVector3d bc) {
        return bc[0]*density_[F_(fi, 0)] + bc[1]*density_[F_(fi, 1)] + bc[2]*density_[F_(fi, 2)];
    }

    // density over a convex polygon which is the overlap part of the triangle and voronoi cell
    double density_integral(Eigen::MatrixX2d& poly, Eigen::VectorXd d);

private:
    Eigen::MatrixXd UV_;
    Eigen::MatrixX2d samples_;
    Eigen::VectorXd density_;
    const Eigen::MatrixXd &V_;
    const Eigen::MatrixX3i &F_;
};
