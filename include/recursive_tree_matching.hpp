// Copyright (c) 2026 Marek Galazka. SPDX-License-Identifier: MIT
#ifndef RECURSIVE_TREE_MATCHING_HPP
#define RECURSIVE_TREE_MATCHING_HPP

#include <algorithm>
#include <cassert>
#include <climits>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace recursive_matching {

struct UpdateStats {
    int new_vertex;
    int changed_states;       // K_n: old vertices only.
    int inspected_vertices;   // U_n: the upward propagation loop only.
    int edge_recourse;        // Symmetric difference of the explicit matchings.
    bool root_changed;
};

// Vertices are numbered 0,1,... in arrival order; construction creates root 0.
// Correctness holds for every leaf-insertion sequence. The paper's stochastic
// bounds additionally require independent, uniform parent choices.
class RecursiveTreeMatching {
public:
    explicit RecursiveTreeMatching(std::size_t reserve_vertices = 1) {
        if (reserve_vertices > static_cast<std::size_t>(INT_MAX))
            throw std::length_error("Too many vertices for int identifiers");
        parent_.reserve(reserve_vertices); head_.reserve(reserve_vertices);
        previous_.reserve(reserve_vertices); next_.reserve(reserve_vertices);
        choice_.reserve(reserve_vertices); mate_.reserve(reserve_vertices);
        state_.reserve(reserve_vertices); pending_.reserve(64);
        append_vertex(-1);
    }

    int vertex_count() const { return static_cast<int>(parent_.size()); }
    int matching_size() const { return cardinality_; }

    int mate_of(int vertex) const {
        check_vertex(vertex);
        return mate_[vertex]; // -1 denotes an unmatched vertex.
    }

    bool contains_edge(int u, int v) const {
        check_vertex(u); check_vertex(v);
        return mate_[u] == v;
    }

    int state_of(int vertex) const {
        check_vertex(vertex);
        return state_[vertex]; // State zero does not mean unmatched.
    }

    int parent_of(int vertex) const {
        check_vertex(vertex);
        return parent_[vertex]; // -1 at the root.
    }

    int selected_child_of(int vertex) const {
        check_vertex(vertex);
        return choice_[vertex];
    }

    UpdateStats add_leaf(int parent) {
        check_vertex(parent);
        if (vertex_count() == INT_MAX)
            throw std::length_error("Too many vertices for int identifiers");
        const int x = vertex_count();
        append_vertex(parent);
        link_zero(x, parent);
        pending_.clear();
        int v = parent;
        UpdateStats result{x, 0, 0, 0, false};
        while (true) {
            ++result.inspected_vertices;
            const int old_state = state_[v];
            state_[v] = head_[v] != -1;
            cardinality_ += static_cast<int>(state_[v]) - old_state;
            const int old_child = choice_[v];
            const int new_child = state_[v] == 0 ? -1 :
                (old_child != -1 && state_[old_child] == 0 ? old_child : head_[v]);
            if (old_child != new_child) {
                pending_.push_back({v, old_child, new_child});
                choice_[v] = new_child;
            }
            if (state_[v] == old_state) break;
            ++result.changed_states;
            if (v == 0) { result.root_changed = true; break; }
            const int q = parent_[v];
            if (state_[v] == 0) link_zero(v, q); else unlink_zero(v, q);
            v = q;
        }
        // Commit every deletion before any addition. Temporary choices are not
        // exposed to clients: queries are made between completed add_leaf calls.
        for (const Choice& change : pending_) {
            if (change.old_child != -1) {
                assert(mate_[change.vertex] == change.old_child);
                assert(mate_[change.old_child] == change.vertex);
                mate_[change.vertex] = mate_[change.old_child] = -1;
                ++result.edge_recourse;
            }
        }
        for (const Choice& change : pending_) {
            if (change.new_child != -1) {
                assert(mate_[change.vertex] == -1 && mate_[change.new_child] == -1);
                mate_[change.vertex] = change.new_child;
                mate_[change.new_child] = change.vertex;
                ++result.edge_recourse;
            }
        }
        assert(result.inspected_vertices == result.changed_states + 1 - result.root_changed);
        assert(result.edge_recourse <= result.changed_states + 2);
        return result;
    }

    // O(n) diagnostic, outside the timed propagation loop. Compute the optimum
    // independently with the static free-root / unrestricted-root matching DP.
    // Assertions must be enabled to use this check.
    void verify_invariants() const {
        const int n = vertex_count();
        std::vector<int> free_value(n, 0), best(n, 0), gain(n, 0), seen(n, 0);
        int edges = 0;
        for (int v = n - 1; v >= 0; --v) {
            best[v] = free_value[v] + gain[v];
            assert(gain[v] == state_[v]);
            if (v > 0) {
                const int p = parent_[v];
                assert(0 <= p && p < v);
                free_value[p] += best[v];
                gain[p] = std::max(gain[p], 1 + free_value[v] - best[v]);
            }
            assert((choice_[v] != -1) == (state_[v] != 0));
            if (choice_[v] != -1) {
                const int c = choice_[v];
                assert(parent_[c] == v && state_[c] == 0 && mate_[v] == c);
                ++edges;
            }
            if (mate_[v] != -1) {
                const int m = mate_[v];
                assert(mate_[m] == v && (parent_[v] == m || parent_[m] == v));
            }
            const int expected_mate = state_[v] ? choice_[v] :
                (v > 0 && choice_[parent_[v]] == v ? parent_[v] : -1);
            assert(expected_mate == mate_[v]);
        }
        assert(edges == cardinality_ && cardinality_ == best[0]);
        for (int p = 0; p < n; ++p) {
            int previous_node = -1;
            for (int v = head_[p]; v != -1; v = next_[v]) {
                assert(parent_[v] == p && state_[v] == 0 && seen[v] == 0);
                assert(previous_[v] == previous_node);
                seen[v] = 1;
                previous_node = v;
            }
        }
        for (int v = 1; v < n; ++v) assert(seen[v] == (state_[v] == 0));
    }

private:
    struct Choice { int vertex, old_child, new_child; };
    std::vector<int> parent_, head_, previous_, next_, choice_, mate_;
    std::vector<unsigned char> state_;
    std::vector<Choice> pending_;
    int cardinality_ = 0;

    void check_vertex(int vertex) const {
        if (vertex < 0 || vertex >= vertex_count())
            throw std::out_of_range("Vertex does not exist");
    }

    void append_vertex(int parent) {
        parent_.push_back(parent); head_.push_back(-1); previous_.push_back(-1);
        next_.push_back(-1); choice_.push_back(-1); mate_.push_back(-1);
        state_.push_back(0);
    }

    void link_zero(int vertex, int parent) {
        assert(state_[vertex] == 0 && previous_[vertex] == -1 && next_[vertex] == -1);
        const int old_head = head_[parent];
        next_[vertex] = old_head;
        if (old_head != -1) previous_[old_head] = vertex;
        head_[parent] = vertex;
    }

    void unlink_zero(int vertex, int parent) {
        if (previous_[vertex] == -1) {
            assert(head_[parent] == vertex);
            head_[parent] = next_[vertex];
        } else next_[previous_[vertex]] = next_[vertex];
        if (next_[vertex] != -1) previous_[next_[vertex]] = previous_[vertex];
        previous_[vertex] = next_[vertex] = -1;
    }
};

} // namespace recursive_matching
#endif
