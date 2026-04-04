#pragma once
#include <vector>
#ifndef ELEMENTS
#define ELEMENTS
    struct Node{
        int id;
        std::vector<double> coords;
    };

    struct Element{
        int id;
        int node_i;
        int node_j;
    };

    struct Support{
        int node_id;
        std::vector<bool> constrained_dofs;
    };

    struct Load{
        int node_id;
        std::vector<double> force;
    };
struct ElementResponse{
    double force;
    double strain;
    double stress;
};
#endif