// Reproducible checks and finite-size experiments for recursive_matching.tex.
// C++17; no external dependencies. Run without NDEBUG so checks stay enabled.
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "recursive_tree_matching.hpp"
using Matching = recursive_matching::RecursiveTreeMatching;
using Update = recursive_matching::UpdateStats;

// mt19937_64 is specified by the C++ standard. Explicit rejection sampling
// avoids implementation-dependent uniform_int_distribution sequences.
std::uint64_t uniform_parent(std::mt19937_64& engine, std::uint64_t bound) {
    const std::uint64_t threshold = (-bound) % bound;
    std::uint64_t sample;
    do sample = engine(); while (sample < threshold);
    return sample % bound;
}

struct Exact {
    std::uint64_t trees = 0, k = 0, u = 0, power_two = 0, root = 0, recourse = 0;
};

void enumerate(int vertex, std::vector<int>& parents, Exact& result) {
    if (vertex == static_cast<int>(parents.size())) {
        Matching matching(static_cast<int>(parents.size()));
        Update last{-1, 0, 0, 0, false};
        for (int x = 1; x < static_cast<int>(parents.size()); ++x) {
            last = matching.add_leaf(parents[x]);
            matching.verify_invariants();
        }
        ++result.trees;
        result.k += last.changed_states; result.u += last.inspected_vertices;
        result.power_two += std::uint64_t(1) << last.changed_states;
        result.root += last.root_changed; result.recourse += last.edge_recourse;
        return;
    }
    for (int p = 0; p < vertex; ++p) {
        parents[vertex] = p;
        enumerate(vertex + 1, parents, result);
    }
}

int main(int argc, char** argv) {
    const std::string directory = argc > 1 ? argv[1] : ".";
    std::ofstream exact_file(directory + "/matching_exact_checks.csv");
    std::ofstream runs_file(directory + "/matching_runs.csv");
    std::ofstream histogram_file(directory + "/matching_late_histogram.csv");
    if (!exact_file || !runs_file || !histogram_file) {
        std::cerr << "Cannot create output files in " << directory << '\n';
        return 1;
    }
    exact_file << "old_n,trees,sum_K,sum_U,sum_2_power_K,root_flips,sum_recourse\n";
    // The before/after example in Fig. 1, including its selected edges.
    Matching illustrated(6);
    illustrated.add_leaf(0); illustrated.add_leaf(0);
    illustrated.add_leaf(2); illustrated.add_leaf(3);
    assert(illustrated.selected_child_of(0) == 1 && illustrated.selected_child_of(3) == 4);
    const Update illustrated_update = illustrated.add_leaf(4);
    illustrated.verify_invariants();
    assert(illustrated_update.changed_states == 3 && illustrated_update.inspected_vertices == 4);
    assert(illustrated_update.edge_recourse == 3 && !illustrated_update.root_changed);
    assert(illustrated.selected_child_of(0) == 1 && illustrated.selected_child_of(2) == 3 && illustrated.selected_child_of(4) == 5);
    std::uint64_t factorial = 1, checked = 1;
    for (int n = 1; n <= 8; ++n) {
        factorial *= n;
        std::vector<int> parents(n + 1, -1);
        Exact exact;
        enumerate(1, parents, exact);
        assert(exact.trees == factorial);
        assert(exact.k == exact.trees);
        assert(n * exact.u == (2 * n - 1) * exact.trees);
        assert(n * exact.power_two == (3 * n - 1) * exact.trees);
        assert(n * exact.root == exact.trees);
        assert(exact.recourse <= 3 * exact.trees);
        exact_file << n << ',' << exact.trees << ',' << exact.k << ',' << exact.u
                   << ',' << exact.power_two << ',' << exact.root << ',' << exact.recourse << '\n';
        checked += exact.trees;
    }
    std::cout << "Full matching implementation: " << checked
              << " recursive trees (sizes 1..9), all checks passed.\n" << std::flush;

    runs_file << "seed,N,S_over_N,W,W_over_ell,recourse_per_update,late_updates,late_mean_K,late_variance_K,late_P_K_zero\n";
    runs_file << std::setprecision(17);
    const std::vector<int> checkpoints{1000, 10000, 100000, 1000000};
    const int final_n = checkpoints.back(), runs = 32;
    std::vector<std::uint64_t> pooled_histogram(64, 0);
    for (int seed = 0; seed < runs; ++seed) {
        std::mt19937_64 engine(static_cast<std::uint64_t>(seed));
        Matching matching(final_n);
        std::uint64_t total_u = 0, total_recourse = 0, late_count = 0, late_zero = 0;
        long double late_sum = 0, late_square_sum = 0;
        int maximum_u = 0;
        std::size_t checkpoint = 0;
        for (int x = 1; x < final_n; ++x) {
            const Update update = matching.add_leaf(static_cast<int>(uniform_parent(engine, x)));
            total_u += update.inspected_vertices; total_recourse += update.edge_recourse;
            maximum_u = std::max(maximum_u, update.inspected_vertices);
            if (x >= final_n / 2) {
                ++late_count; late_sum += update.changed_states;
                late_zero += update.changed_states == 0;
                late_square_sum += update.changed_states * update.changed_states;
                assert(update.changed_states < static_cast<int>(pooled_histogram.size()));
                ++pooled_histogram[update.changed_states];
            }
            if (x + 1 == checkpoints[checkpoint]) {
                matching.verify_invariants();
                const int n = x + 1;
                const double ell = std::log(static_cast<double>(n)) / std::log(std::log(static_cast<double>(n)));
                runs_file << seed << ',' << n << ',' << static_cast<double>(total_u) / n
                          << ',' << maximum_u << ',' << maximum_u / ell << ','
                          << static_cast<double>(total_recourse) / (n - 1) << ',' << late_count;
                if (late_count) {
                    const long double mean = late_sum / late_count;
                    runs_file << ',' << static_cast<double>(mean) << ','
                              << static_cast<double>(late_square_sum / late_count - mean * mean)
                              << ',' << static_cast<double>(late_zero) / late_count;
                } else runs_file << ",,,";
                runs_file << '\n';
                ++checkpoint;
            }
        }
        std::cout << "seed " << seed << " finished, W=" << maximum_u << '\n' << std::flush;
    }
    histogram_file << "K,count\n";
    for (std::size_t k = 0; k < pooled_histogram.size(); ++k)
        if (pooled_histogram[k]) histogram_file << k << ',' << pooled_histogram[k] << '\n';
    std::cout << "All simulations and checkpoint matching verifications passed.\n";
}
