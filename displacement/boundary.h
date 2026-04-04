#pragma once
#include "../components/element.h"
#include "../components/solver.h"
#include <vector>

void inline applyBoundaryConditions(
    Matrix& K_global,
    std::vector<double>& F,
    const std::vector<Node>& nodes,
    const std::vector<Support>& supports,
    const int dim   
){
    const int n_nodes = nodes.size();
    const int n_dofs = n_nodes*dim;
    for(const auto& sup: supports){
    int offset = sup.node_id * dim;
    for(int d=0;d<dim;d++){
        if(sup.constrained_dofs[d]){
            int fixed_dof = offset + d;
            for(int j=0;j<n_dofs;j++){
                K_global(fixed_dof,j) = 0.0;
                K_global(j,fixed_dof) = 0.0; 
            }
            K_global(fixed_dof,fixed_dof)=1.0;
            F[fixed_dof] = 0.0;
        }
    }
}
}