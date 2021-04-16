//
// Created by pupa on 4/16/21.
//

#pragma once


#include <complex>
#include <vector>

using Complex = std::complex<double>;

using namespace std::complex_literals;

bool inline to_left_test(const Complex& p, const Complex& p1, const Complex& p2) {
    return std::imag(std::conj(p2-p1)*(p-p1)) > 0;
}

Complex inline intersection(const Complex& p0, Complex v0, const Complex& p1, Complex v1) {
    Complex _n0 = std::conj(v0 * 1.0i /std::abs(v0));
    // dot(n0, p1+t1*v1) = dot(n0, p0)
    double _t1 = std::real(_n0*p0 - _n0*p1)/std::real(_n0*v1);
    return p1 + _t1 * v1;
}

// Sutherland-Hodgman clipping
void static clipping_sutherland_hodgman(std::vector<Complex>& subj, std::vector<Complex>& clip, std::vector<Complex>& solution) {
    Complex cp1, cp2, s, e;

    std::vector<Complex> input;
    solution = subj;

    for(int j = 0; j < clip.size(); j++) {
        std::swap(input, solution);
        solution.clear();

        // get clipping polygon edge
        cp1 = clip[j]; cp2 = clip[(j + 1) % clip.size()];

        for(int i = 0; i < input.size(); i++) {
            // get subject polygon edge
            s = input[i]; e = input[(i + 1) % input.size()];

            // Case 1: Both vertices are inside:
            // Only the second vertex is added to the output list
            if(to_left_test(s, cp1, cp2) && to_left_test(e, cp1, cp2))
                solution.push_back(e);
                // Case 2: First vertex is outside while second one is inside:
                // Both the point of intersection of the edge with the clip boundary
                // and the second vertex are added to the output list
            else if(!to_left_test(s, cp1, cp2) && to_left_test(e, cp1, cp2)){
                solution.push_back( intersection(cp1, cp2-cp1, s, e-s) );
                solution.push_back(e);
            }
                // Case 3: First vertex is inside while second one is outside:
                // Only the point of intersection of the edge with the clip boundary
                // is added to the output list
            else if(to_left_test(s, cp1, cp2) && !to_left_test(e, cp1, cp2))
                solution.push_back(  intersection(cp1, cp2 - cp1, s, e -s) );
                // Case 4: Both vertices are outside
                // No vertices are added to the output list
            else if(!to_left_test(s, cp1, cp2) && !to_left_test(e, cp1, cp2))
                ;
        }

    }
}