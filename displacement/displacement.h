#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <utility>
#include <numeric>
#include "../components/solver.h"
#include "../components/element.h"
#include "assembly.h"
#include "boundary.h"

std::vector<double> solve_disp_truss(
    const std::vector<Node>& nodes,
    const std::vector<Element>& elements,
    const std::vector<Support>& supports,
    const std::vector<Load>& loads,
    int dim,
    const double E,
    const double A
){
    int n_nodes = nodes.size();
    int n_dofs = n_nodes*dim;

    Matrix K_global(n_dofs,n_dofs);
    std::vector<double> F(n_dofs,0.0);

    for(const auto& load:loads){
        int base = load.node_id * dim;
        for(int d=0;d<dim;d++){
            F[base + d] += load.force[d];
        }
    }

    double stiff_coeff = E*A;
    assembleGlobalStiffness(K_global,nodes,elements,dim,stiff_coeff);

    applyBoundaryConditions(K_global,F,nodes,supports,dim);
    for(int i=0;i<n_dofs;i++){
    if(std::abs(K_global(i,i)) < 1e-12){
        std::cout << "Near-zero diagonal at DOF " << i << "\n";
    }
}
    std::vector<double> u;
    conjugate_gradient(K_global,F,u);

    return u;
}