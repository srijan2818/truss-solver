#pragma once
#include "../components/element.h"
#include "../components/solver.h"
#include <vector>

Matrix elementStiffnessMatrix(
    const Node& node1, 
    const Node& node2,
    const int dim,
    double stiff_coeff
){
    double length = 0;

    double n[3]; // it would be most likely at most 3, access first 'dim' elements, reduces O(num_elements) heap ops
    for(int d=0;d<dim;d++){
                n[d] = node2.coords[d] - node1.coords[d]; 
                length += n[d] * n[d];
            } 
    length = std::sqrt(length);  
    if(length < 1e-12){
    throw std::runtime_error("Zero-length element detected");
}
    for(int d=0;d<dim;d++)
    {
        n[d] /= length;
    }
    stiff_coeff /= length;

    Matrix K_local(2*dim,2*dim);
    for(int row=0;row<2;row++){
        for(int col=0;col<2;col++){
            double sign = (row==col) ? 1.0 : -1.0;
            for(int r=0;r<dim;r++){
                for(int c=0;c<dim;c++){
                    int local_i = row*dim +r;
                    int local_j = col*dim + c;
                    K_local(local_i,local_j) += stiff_coeff * n[r] * n[c] * sign;
                }
        }
    }
    }
    return K_local;
}

inline void assembleGlobalStiffness(
    Matrix& K_global,
    const std::vector<Node>& nodes,
    const std::vector<Element>& elements,
    const int dim,
    double stiff_coeff
){
    const int n_nodes = nodes.size();
    const int n_dofs = n_nodes*dim;
    
    for(const auto& element:elements){
        const int node_i = element.node_i;
        const int node_j = element.node_j;

        const int node_ids[2] = {node_i,node_j};
    double length = 0;

    std::vector<double> n(3);
    for(int d=0;d<dim;d++){
                n[d] = nodes[node_j].coords[d] - nodes[node_i].coords[d]; 
                length += n[d] * n[d];
            } 
    length = std::sqrt(length);
    if(length < 1e-12){
    throw std::runtime_error("Zero-length element detected");
}
    for(int d=0;d<dim;d++)
    {
        n[d] /= length;
    }
    double k = stiff_coeff / length;
    for(int r=0; r<dim;r++){
        int offset_i = node_i*dim;
        int offset_j = node_j*dim;
        int I_i = offset_i+ r;
        int I_j = offset_j + r;

        for(int c=0;c<dim;c++){
            int J_i = offset_i + c;
            int J_j = offset_j + c;
            double val = k * n[r] * n[c];
            K_global(I_i, J_i) +=  val;  // +nnT
            K_global(I_i, J_j) += -val;  // -nnT
            K_global(I_j, J_i) += -val;  // -nnT
            K_global(I_j, J_j) +=  val;  // +nnT
        }
    }     
    }

}
inline void elementDofMap(
    int* map,
    const int node_i,
    const int node_j,
    const int dim
){
    for(int d=0;d<dim;d++){
        map[d] = node_i*dim + d;
        map[dim + d] = node_j*dim + d;
    }
}