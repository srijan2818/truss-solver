#include <vector>
#include <iostream>
#include <string>
#include <unordered_set>
#include <stack>
#include "../components/element.h"

struct ValidationResult{
    bool result;
    bool is_determinate;
    std::string msg;
};

// 1. Length checks inside solvers to avoid recomputation
// 2. Duplicate element check
inline long long encode(int u,int v){
    return (static_cast<long long>(u) << 32) | v;
}
inline bool checkDuplicateElements(
    const std::vector<Element>& elements
){
    const int n_elements = elements.size();
    std::unordered_set<long long> seen;
    seen.reserve(n_elements);
    for(const auto& elem:elements){
        const int i = elem.node_i;
        const int j = elem.node_j;

        int u = std::min(i,j);
        int v = std::max(i,j);

        long long key = encode(u,v);

        if(seen.find(key) != seen.end()){
            std::cout << "Duplicate element found between " << u << " - " << v << "\n";
            return false;
        }

        seen.insert(key);
    }
    return true;

}
// 3. Inalide node index
inline bool checkElementNodeValidity(
    const std::vector<Node>& nodes,
    const std::vector<Element>& elements
){
    const int n = nodes.size();

    for(const auto& elem : elements){
        const int i = elem.node_i;
        const int j = elem.node_j;
        if(i < 0 || i >= n || j < 0 || j >= n){
            std::cout << "Invalid node index in element\n";
            return false;
        }
    }

    return true;
}
// 4. TOPOLOGY FAILS - start from supports nd check connectivity everywhere
inline bool connectivityCheck(
    const std::vector<Node>&nodes,
    const std::vector<Element>& elements,
    const int dim,
    const std::vector<Support>& supports
){
    int n_nodes = nodes.size();
    int n_elements = elements.size();

    if (n_nodes == 0 || n_elements ==0 || supports.empty()) return false;

    std::vector<int> recent_edge_id (n_nodes,-1);
    std::vector<int> to   (2*n_elements);
    std::vector<int> last_edge_id (2*n_elements);
    int edgePtr = 0;
    for(const auto& elem:elements){

        int i = elem.node_i;
        int j = elem.node_j;

        to[edgePtr] = j;
        last_edge_id[edgePtr] = recent_edge_id[i];
        recent_edge_id[i] = edgePtr++;

        to[edgePtr] = i;
        last_edge_id[edgePtr] = recent_edge_id[j];
        recent_edge_id[j] = edgePtr++;
        
    }

    std::vector<bool> visited(n_nodes,false);
    std::stack<int> s;
    for (const auto& sup : supports) {
        if (!visited[sup.node_id]) {
            s.push(sup.node_id);
            visited[sup.node_id] = true;
        }
    }
    
    if (s.empty()) return false;
    int countVisited = 0;

    while(!s.empty()){
        int u = s.top();
        s.pop();
        countVisited++;

        for(int e = recent_edge_id[u]; e!= -1; e = last_edge_id[e]){
            int v = to[e];
            if(!visited[v]){
                visited[v] = true;
                s.push(v);
            }
        }
    }
    return countVisited == n_nodes;

}

// 5. Constraint validation - dof nd solver lvl
inline int dof_constraint(
    const std::vector<Node>&nodes,
    const std::vector<Element>& elements,
    const int dim,
    const std::vector<Support>& supports
){
        int n_nodes = nodes.size();
        int n_elements = elements.size();
        int n_eqns = n_nodes * dim;

        int n_unknowns = n_elements;
        for(const auto& sup:supports){
            for(bool constrained: sup.constrained_dofs){
                if(constrained) n_unknowns++;
            }
        }
        if(n_eqns < n_unknowns) return -1;        // Statically Indeterminate - overconstrained
        else if (n_eqns == n_unknowns) return 0;  // Statically Determinate
        else return 1;                            // mechanism - unstable/ underconstrained
}


// 6. Material check
inline bool checkMaterial(const double E,const double A){
    if(E<=0 || A<=0){
        std::cout<<"Invalid Material properties  ( E or A <= 0 )\n";
        return false;
    }
    return true;
}


// 7. Loads check
inline bool checkLoads(
    const std::vector<Node>& nodes,
    const std::vector<Load>& loads,
    int dim
){
    int n = nodes.size();

    for(const auto& l : loads){
        if(l.node_id < 0 || l.node_id >= n){
            std::cout << "Invalid load node index\n";
            return false;
        }

        if(l.force.size() != dim){
            std::cout << "Load dimension mismatch\n";
            return false;
        }
    }

    return true;
}

inline ValidationResult verifyModel(
    const std::vector<Node>& nodes,
    const std::vector<Element>& elements,
    const std::vector<Support>& supports,
    const std::vector<Load>& loads,
    const int dim,
    const double E,
    const double A
){
    if(!checkElementNodeValidity(nodes,elements)){
        return {false,false,"Invalid node indices"};
    }
    if(!checkDuplicateElements(elements)){
        return {false,false,"Duplicate elements detected"};
    }
    if(!connectivityCheck(nodes,elements,dim,supports)){
        return {false,false,"Disconnected structure"};
    }
    if(!checkMaterial(E,A)){
        return {false,false,"Invalid material properties"};
    }
    if(!checkLoads(nodes,loads,dim)){
        return {false,false,"Invalid loads"};
    }
    int type = dof_constraint(nodes,elements,dim,supports);

    if(type == 1){
        return {false,false,"Mechanism (under-constrained system)"};
    }

    bool is_determinate = (type == 0);

    return {true, is_determinate, "OK"};

}