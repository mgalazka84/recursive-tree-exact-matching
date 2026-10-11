// The insertion shown in Figure 1. Code labels are paper labels minus one.
#include "recursive_tree_matching.hpp"
#include <iostream>

int main() {
    recursive_matching::RecursiveTreeMatching tree(6);
    tree.add_leaf(0); tree.add_leaf(0); tree.add_leaf(2); tree.add_leaf(3);
    const auto update = tree.add_leaf(4);
    tree.verify_invariants();
    std::cout << "K=" << update.changed_states << ", U=" << update.inspected_vertices
              << ", recourse=" << update.edge_recourse << '\n';
    std::cout << "Maximum matching size: " << tree.matching_size() << '\n';
    for (int v = 0; v < tree.vertex_count(); ++v) {
        std::cout << "vertex " << v << ": state=" << tree.state_of(v)
                  << ", mate=" << tree.mate_of(v) << '\n';
    }
    std::cout << "Edge (0,1) is selected: " << tree.contains_edge(0,1) << '\n';
}
