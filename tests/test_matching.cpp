// Independent edge-subset enumeration, plus the public query contracts.
#include "recursive_tree_matching.hpp"
#include <cassert>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>

int brute_force_optimum(const std::vector<int>& parents, int n) {
    int optimum = 0;
    for (unsigned subset = 0; subset < (1u << (n - 1)); ++subset) {
        std::vector<bool> used(n, false);
        int count = 0;
        bool valid = true;
        for (int v = 1; v < n; ++v) if (subset & (1u << (v - 1))) {
            const int p = parents[v];
            if (used[v] || used[p]) { valid = false; break; }
            used[v] = used[p] = true;
            ++count;
        }
        if (valid) optimum = std::max(optimum, count);
    }
    return optimum;
}

std::uint64_t checked = 0;
void enumerate(int v, std::vector<int>& parents) {
    if (v == static_cast<int>(parents.size())) {
        recursive_matching::RecursiveTreeMatching tree;
        assert(tree.vertex_count() == 1 && tree.matching_size() == 0);
        for (int n = 1; n <= static_cast<int>(parents.size()); ++n) {
            if (n > 1) {
                const auto result = tree.add_leaf(parents[n - 1]);
                assert(result.new_vertex == n - 1);
            }
            tree.verify_invariants();
            assert(tree.matching_size() == brute_force_optimum(parents, n));
            int endpoints = 0;
            for (int a = 0; a < n; ++a) {
                const int partner = tree.mate_of(a);
                endpoints += partner != -1;
                if (partner != -1) {
                    assert(tree.mate_of(partner) == a);
                    assert(parents[a] == partner || parents[partner] == a);
                }
                for (int b = 0; b < n; ++b)
                    assert(tree.contains_edge(a,b) == (partner == b));
            }
            assert(endpoints == 2 * tree.matching_size());
        }
        ++checked;
        return;
    }
    for (int p = 0; p < v; ++p) {
        parents[v] = p;
        enumerate(v + 1, parents);
    }
}

int main() {
    for (int n = 1; n <= 7; ++n) {
        std::vector<int> parents(n, -1);
        enumerate(1, parents);
    }
    recursive_matching::RecursiveTreeMatching one;
    bool invalid_parent = false, invalid_query = false;
    try { one.add_leaf(1); } catch (const std::out_of_range&) { invalid_parent = true; }
    try { one.mate_of(-1); } catch (const std::out_of_range&) { invalid_query = true; }
    assert(invalid_parent && invalid_query && one.vertex_count() == 1);
    assert(one.parent_of(0) == -1 && !one.contains_edge(0,0));
    std::cout << "All public queries and brute-force optima passed on "
              << checked << " recursive trees of sizes 1..7.\n";
}
