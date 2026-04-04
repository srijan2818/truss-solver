#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <utility>
#include <numeric>
#include "../components/solver.h"
#include "../components/element.h"
    
    std::vector<double> solve_truss(const std::vector<Node>& nodes, const std::vector<Element>& elements,const std::vector<Support>& supports, const std::vector<Load>& loads, int dim){
        
        int n_nodes = nodes.size();
        int n_elements = elements.size();
        int n_eqns = n_nodes * dim;

        int n_unknowns = n_elements;
        for(const auto& sup:supports){
            for(bool constrained: sup.constrained_dofs){
                if(constrained) n_unknowns++;
            }
        }

        if (n_unknowns != n_eqns){
            std::cout << "System is not statically determinable\n Number of unknowns : " << n_unknowns << "\nNumber of equations : \n" << n_eqns;
            throw std::runtime_error("System not statically determinate");
        }

        Matrix A(n_eqns,n_unknowns);
        std::vector<double> b(n_eqns,0.0);
        for(int k=0;k<n_elements;k++){
            const auto& element = elements[k];
            int n_i = element.node_i;
            int n_j = element.node_j;

            double length = 0.0;
            std::vector<double> dir_vec(dim);
            for (int d=0;d<dim;d++){
                dir_vec[d] = nodes[n_j].coords[d] - nodes[n_i].coords[d];
                length += dir_vec[d] * dir_vec[d];
            }

            length = std::sqrt(length);
            if (length < 1e-12)
               throw std::runtime_error("Zero length element");
            for(int d=0;d<dim;d++){
                double c = dir_vec[d]/length;
                A(n_i*dim+d,k) = c;
                A(n_j*dim+d,k) = -c;

            }
        }

        int reaction_col = n_elements;
        for(const auto& sup:supports){
            for(int d=0;d<dim;d++){
                if(sup.constrained_dofs[d]){
                    A(sup.node_id*dim + d,reaction_col) = 1.0;
                    reaction_col++;
                }
            }
        }

        for(const auto& load:loads){
            for(int d=0;d<dim;d++){
                b[load.node_id*dim +d] += -load.force[d];
            }
        }


        //Verify
    std::vector<double> solutions;
    // lu_solve(A, b, solutions, "doolittle");
    qr_solve(A,b,solutions,"householder");
    bool verified = true;
    for(int i =0;i<n_eqns;i++){
        double res = 0.0;
        for(int j=0;j<n_unknowns;j++){
            res += A(i,j) * solutions[j];
        }
        if(std::abs(res-b[i])>1e-6){
            verified=false;
            std::cout<<"Eqll violated at node"<<i/dim<<" residual of "<<std::abs(res-b[i]);
        }
    }
    if(!verified) throw std::runtime_error("Verification failure");
    else std::cout<<"VERIFIED\n";
    return solutions;
    }