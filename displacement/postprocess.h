#pragma once

#include <vector>
#include <cmath>
#include "../components/element.h"

inline void computeElementKinematics(
    const Node& node_i,
    const Node& node_j,
    const std::vector<double>& u,
    const int dim,
    double& delta,
    double& length
){
    double n[3];
    length = 0.0;

    for(int d=0;d<dim;d++){
        n[d] = node_j.coords[d] - node_i.coords[d];
        length += n[d]*n[d];
    }

    length = std::sqrt(length);
    if(length < 1e-12){
    throw std::runtime_error("Zero-length element detected");
}

    for(int d=0;d<dim;d++){
        n[d] /= length;
    }
    const int idx_i = node_i.id;
    const int idx_j = node_j.id; 
    int base_i = idx_i * dim;
    int base_j = idx_j * dim;

    delta = 0.0;
    for(int d=0;d<dim;d++){
        delta += n[d] * (u[base_j + d] - u[base_i + d]);
    }
}

inline double computeElementForce(
    const Node& node_i,
    const Node& node_j,
    const std::vector<double>& u,
    const int dim,
    const double stiff_coeff
){
    double delta,length;
    computeElementKinematics(node_i,node_j,u,dim,delta,length);

    double f = (stiff_coeff/length) * delta;

    return f;
}

inline double computeElementStrain(
    const Node& node_i,
    const Node& node_j,
    const std::vector<double>& u,
    const int dim
){
    double delta,length;
    computeElementKinematics(node_i,node_j,u,dim,delta,length);

    return delta/length;
}

inline double computeElementStress(
    const Node& node_i,
    const Node& node_j,

    const std::vector<double>& u,
    const int dim,
    const double E
){
    double strain = computeElementStrain(node_i,node_j,u,dim);

    return E*strain;
}

inline ElementResponse computeElementResponse(
    const Node& node_i,
    const Node& node_j,
    const std::vector<double>& u,
    const int dim,
    const double E,
    const double A
){
    double delta,length;
    computeElementKinematics(node_i,node_j,u,dim,delta,length);

    double strain = delta/length;
    double stress = E*strain;
    double force  = A*stress;

    return {force,strain,stress};
}

inline std::vector<double> computeAllElementForces(
    const std::vector<Node>& nodes,
    const std::vector<Element>& elements,
    const std::vector<double>& u,
    const int dim,
    const double E,
    const double A
){
    std::vector<double> forces(elements.size());

    for(int e=0;e<elements.size();e++){
        const auto& element = elements[e];

        const Node& node_i = nodes[element.node_i];
        const Node& node_j = nodes[element.node_j];
        const double stiff_coeff = E*A;

        forces[e] = computeElementForce(node_i,node_j,u,dim,stiff_coeff);
    }

    return forces;
}

inline std::vector<ElementResponse> computeAllElementResponses(
    const std::vector<Node>& nodes,
    const std::vector<Element>& elements,
    const std::vector<double>& u,
    const int dim,
    const double E,
    const double A
){
    std::vector<ElementResponse> results(elements.size());

    for(int e=0;e<elements.size();e++){
        const auto& element = elements[e];

        results[e] = computeElementResponse(nodes[element.node_i],nodes[element.node_j],u,dim,E,A);
    }

    return results;
}

