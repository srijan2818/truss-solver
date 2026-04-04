#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>

#include "./components/element.h"
#include "./verification/validate.h"
#include "./displacement/displacement.h"
#include "./displacement/postprocess.h"
#include "./equilibrium/statically_determined.h"

// Simple reporter
void report(std::string name, bool pass, std::string details = "") {
    std::cout << std::setw(35) << std::left << name 
              << " [" << (pass ? "PASS" : "FAIL") << "] " 
              << details << "\n";
}

int main() {
    std::cout << "--- TRUSS SOLVER TEST SUITE ---\n" << std::endl;

    const int dim = 2;
    const double E = 200e9;
    const double A = 0.01;

    // CASE 1: Mechanism (Under-constrained)
    // Expected: Validation should catch the DOF imbalance
    {
        std::vector<Node> nodes = {{0, {0.0, 0.0}}, {1, {1.0, 0.0}}};
        std::vector<Element> elements = {{0, 0, 1}};
        std::vector<Support> supports = {{0, {true, true}}}; // Node 1 free to rotate
        std::vector<Load> loads = {{1, {0.0, -100.0}}};

        ValidationResult res = verifyModel(nodes, elements, supports, loads, dim, E, A);
        report("Mechanism Check", !res.result, "Caught under-constrained system");
    }

    // CASE 2: Disconnected Node
    // Expected: Connectivity check should fail
    {
        std::vector<Node> nodes = {{0, {0.0, 0.0}}, {1, {1.0, 0.0}}, {2, {2.0, 0.0}}};
        std::vector<Element> elements = {{0, 0, 1}}; // Node 2 is floating
        std::vector<Support> supports = {{0, {true, true}}, {1, {false, true}}};
        std::vector<Load> loads = {{1, {100.0, 0.0}}};

        ValidationResult res = verifyModel(nodes, elements, supports, loads, dim, E, A);
        report("Disconnected Node Check", !res.result, "Caught floating node 2");
    }

    // CASE 3: Single Bar Analytical (Direct delta = FL/EA)
    // Expected: UX = 5.0e-5 for 100kN, 1m bar, 200GPa, 0.01 area
    {
        std::vector<Node> nodes = {{0, {0.0, 0.0}}, {1, {1.0, 0.0}}};
        std::vector<Element> elements = {{0, 0, 1}};
        std::vector<Support> supports = {{0, {true, true}}, {1, {false, true}}};
        std::vector<Load> loads = {{1, {100000.0, 0.0}}};

        auto u = solve_disp_truss(nodes, elements, supports, loads, dim, E, A);
        double expected = 100000.0 * 1.0 / (E * A);
        double diff = std::abs(u[2] - expected);
        
        report("Single Bar Analytical", diff < 1e-12, "Error: " + std::to_string(diff));
    }

    // CASE 4: Warren Truss Cross-Verification
    // Expected: Displacement Method Force == Static Equilibrium Force
    {
        std::vector<Node> nodes = {
            {0, {0.0, 0.0}}, {1, {5.0, 5.0}}, {2, {10.0, 0.0}}
        };
        std::vector<Element> elements = {
            {0, 0, 1}, {1, 1, 2}, {2, 0, 2}
        };
        std::vector<Support> supports = {
            {0, {true, true}}, {2, {false, true}}
        };
        std::vector<Load> loads = {
            {1, {0.0, -100000.0}}
        };

        // Solve both ways
        auto u = solve_disp_truss(nodes, elements, supports, loads, dim, E, A);
        auto res_disp = computeAllElementResponses(nodes, elements, u, dim, E, A);
        auto res_stat = solve_truss(nodes, elements, supports, loads, dim);

        double max_diff = 0;
        for(size_t i=0; i<elements.size(); ++i) {
            max_diff = std::max(max_diff, std::abs(res_disp[i].force - res_stat[i]));
        }

        report("Warren Cross-Verification", max_diff < 1e-6, "Max Force Diff: " + std::to_string(max_diff));
    }

    std::cout << "\n--- TEST COMPLETE ---" << std::endl;
    return 0;
}