# Exact matching in growing trees

Implementation and reproducibility code for two separate manuscripts by
**Marek Gałązka**, Adam Mickiewicz University, Poznań.

| Manuscript | Source | Scope |
| --- | --- | --- |
| **Exact Matching in Growing Recursive Trees with Constant Expected Update Cost** | [paper/recursive_matching.tex](paper/recursive_matching.tex) | Uniform recursive attachment; 12 pages including references and all proofs |
| **Exact Matching under Affine Preferential Attachment: Update Laws, Gaussian Fluctuations, and Extremes** | [paper_affine/affine_matching.tex](paper_affine/affine_matching.tex) | Fixed positive reinforcement; 10 main-text pages, 1 references page, 7 appendix pages |

The first source is restored **byte for byte** from pre-affine commit
`fd72ed52a1eea379ccde2fedcd7cfa2c409d2b00`. The original uniform-experiment
snapshot remains `0e4c10c707699f61ca2e83cd28de8ca4683ef03d`.
The second paper has its own definitions, theorem numbering, proofs, and
bibliography. It explicitly cites the uniform companion and needs none of
its probabilistic theorems. The algorithm implementation is shared.
The earlier combined development version remains in Git history at
`d3b3ace968a14ae8739c061d0c7d451aa516f0d1`.

## Reproduction

Requirements: Python 3.8 or later, a C++17 compiler, and no external Python
packages. LaTeX is required only for the manuscripts. The supplied
reproduction was checked with Python 3.12.14 and GNU C++ 13.3.0.

```sh
python3 reproduce.py --checks-only  # all API, uniform, affine and CLT finite checks
python3 reproduce.py                # also regenerate the uniform simulation table
make demo                          # uniform paper's exact insertion example
make paper                         # uniform manuscript
make paper-affine                  # affine manuscript
make papers                        # both manuscripts
```

New outputs go to `results/` and are compared against `data/`. The full
run regenerates 32 uniform trajectories to one million vertices.
The Python uniform checker retains 12 smaller supporting runs to 100,000
vertices. The affine and CLT checkers use exact arithmetic, not Monte Carlo.
Keep C++ assertions enabled; do not use `-DNDEBUG`. Set `CXX=clang++` to
select another compiler. `make checks` retains its original API and uniform
checks; use `reproduce.py --checks-only` for every model.

## Using the algorithm

Include `include/recursive_tree_matching.hpp`. Vertices are numbered **0,1,...** in arrival order; construction creates the single root 0. The paper instead uses labels 1,2,..., so add one to code identifiers when reading the diagram.

```cpp
#include "recursive_tree_matching.hpp"

recursive_matching::RecursiveTreeMatching tree(6); // optional reserve
tree.add_leaf(0); // create vertex 1, attached to root 0
tree.add_leaf(0); // vertex 2
tree.add_leaf(2); // vertex 3
tree.add_leaf(3); // vertex 4
auto cost = tree.add_leaf(4); // vertex 5: Figure 1's insertion

int cardinality = tree.matching_size(); // 3
int partner = tree.mate_of(0);         // 1; -1 means unmatched
bool selected = tree.contains_edge(0, 1); // true
// cost.changed_states == 3, cost.inspected_vertices == 4,
// cost.edge_recourse == 3, cost.root_changed == false
```

See `examples/demo.cpp` for an executable client. The methods `matching_size()`, `mate_of(v)`, and `contains_edge(u,v)` are constant-time queries between completed insertions. `state_of`, `parent_of`, and `selected_child_of` expose diagnostic information. Invalid vertex identifiers throw `std::out_of_range`.

`add_leaf(parent)` returns `UpdateStats`, including the new identifier and the three counters. The maintained optimum uses zero-child lists, retains an existing valid selected child, and commits all obsolete-edge deletions before new additions. **State zero does not mean unmatched**: a zero-state vertex can be matched to its parent.

Correctness holds for every leaf-insertion history. Each insertion takes
O(1+U_n) propagation time and changes at most K_n+2 matching edges.
Queries take constant time and storage is linear. The parent is supplied;
sampling a random parent is outside the update cost. Edge weights and
deletions are outside the guarantees. The reference implementation uses
`std::vector`; allocation without advance reservation is amortized.
The experiments measure inspections and recourse, not wall-clock time.

## Uniform manuscript

Theorem 1 gives `E K_n = 1`, `E U_n = 2-1/n`, expected recourse at most
three, and uniform exponential-moment bounds. Lemmas 2--3 describe the
known deterministic algorithm; Lemmas 4--6 give sensitivity and fringe
identities. Theorem 7 identifies the limiting update law. Lemma 8 supplies
the fixed-size variance estimate used in Theorem 9's quantitative L1 law.
Theorem 10 proves `W_N/(log N/log log N) -> 1` in probability.

The original figure, simulation table, numbering, and all proofs are
unchanged. This paper does not claim a cumulative-cost CLT. Uniform data
have seeds 0--31, final size 1,000,000, and checkpoints at 10^3, 10^4,
10^5 and 10^6 on each trajectory. Table 1's SD is across independent
trajectories; checkpoints and updates within a trajectory are dependent.
At N=10^6, mean S/N is 1.9998017, mean recourse is 1.307612, and mean
W/(log N/log log N) is 2.126309, with W between 10 and 13. These finite
data do not establish a convergence rate or improve the recourse bound.

## Affine manuscript

Parent weights are `1 + beta * number_of_children`, for fixed `beta > 0`.
The seed-edge convention uses `degree + delta`, `delta > -1`, corresponding
to `beta=1/(1+delta)` in each of its two affine components.

| Result in the affine paper | Statement | Proof |
| --- | --- | --- |
| Lemma 1 | Recalled matching recurrence, implementation and recourse | Section 2 |
| Theorem 2 | Constant expected update cost and factorial update tails | Appendix A |
| Theorem 3 | Parameter-dependent update limit, moments and sharp limiting tails | Appendix B |
| Theorem 4 | `W_N = r_beta(N) + O_beta(1)` almost surely; Lambert-W centering and second-order expansion | Appendix C |
| Theorem 5 | `E S_N = mu_beta N - c_beta + o(1)` with explicit correction | Appendix D |
| Theorem 6 | Nondegenerate cumulative-cost CLT; `Var(S_N) = sigma_beta^2 N + o(N)`; exact N^(-1/2) mean-error order | Section 4 |
| Corollary 7 | Degree-model transfer; ordinary degree attachment has `E U_n <= 3/2` | Appendix E and Section 4 |
| Lemma 8 | Fringe means, cited size-count variance, and independent interiors | Section 3 |
| Lemma 9 | Variance bound for tolls centered at every tree size | Section 4 |

The CLT proof preserves the expected history score at each fringe size.
The residual satisfies `E(A_N-A_tilde_N^(L))^2 <= C_beta N L^(-zeta)`,
where `zeta=1/(1+beta)`. The approximation meets the conditions of the
Ralaivaosaona--Wagner theorem. Four-vertex resampling preserves every
exterior state history and proves the strictly positive lower bound

`sigma_beta^2 >= 2(1+2 beta)(1+beta)^2 / ((2+beta)(3+2 beta)^2(4+3 beta)(5+4 beta))`.

No closed form for sigma_beta^2, beta-uniform bound, Gumbel law, or CLT
at beta=0 is claimed. The weaker concentration estimate from the combined
development draft is omitted because the CLT gives the exact mean-error
order. The affine paper does not reuse the uniform simulation table.

## Code and data

| File | Role and manuscript connection |
| --- | --- |
| `include/recursive_tree_matching.hpp` | Shared complete C++17 matching implementation with size, mate and edge queries |
| `examples/demo.cpp` | Uniform paper, Figure 1 |
| `tests/test_matching.cpp` | Independent edge-subset optimum and query checks on 874 trees of sizes 1--7 |
| `verification/matching_experiments.cpp` | Exact checks on 46,234 trees and the 32 uniform growth trajectories |
| `verification/verify_extensions.py` | Uniform paper's history identities, fringe covariance, finite PGFs and path lower bounds |
| `verification/verify_preferential.py` | Affine paper's sensitivities, update laws, fringe/witness variances, means and component comparisons |
| `verification/verify_clt.py` | Affine Theorem 6 and Lemma 9: centered variance identities and four-vertex replacements |
| `verification/summarize_experiments.py` | Uniform paper's Table 1 |
| `reproduce.py` | Shared one-command reproduction with saved-data comparisons |

Run the affine finite checkers separately with:

```sh
python3 verification/verify_preferential.py --output results/preferential_exact_checks.json
python3 verification/verify_clt.py --output results/clt_exact_checks.json
```

They enumerate all 5,914 trees of sizes 1--8 for seven beta values
`0, 1/10, 1/4, 1/2, 1, 2, 10` and four degree offsets `-1/2, 0, 1, 3`.
The uniform boundary is retained as a finite-code consistency check.
All 46,233 next insertions are checked independently of the sensitivity
recurrence, and formal-series checks run through degree 24. The CLT checker
verifies 252 centered variance identities and 3,474 replacements,
checking every exterior state at every arrival prefix. These finite checks
are not formal proof verification and do not establish asymptotic normality.

Reference outputs are in `data/preferential_exact_checks.json` and
`data/clt_exact_checks.json`. Certified limiting-law intervals use the
proved remainder `zeta^2/(L+zeta)`. The uniform raw files are
`matching_exact_checks.csv`, `matching_runs.csv`, and
`matching_late_histogram.csv`, with `experiment_summary.json` and
`experiment_table.tex`. `verification/verification_results.json` records
the original Python checks and smaller supporting simulations.
Regenerate the uniform table without rerunning trajectories using:

```sh
python3 verification/summarize_experiments.py data
```

## Building and citing

Both manuscripts embed their bibliography and use standard LNCS dimensions.
Their local README files give theorem maps and Overleaf instructions.
The uniform source contains a fallback copy of Table 1. The affine source
needs no external table or figure files. `llncs.cls` retains Springer's
notices and CC BY 4.0 terms; see `THIRD_PARTY_NOTICES.md`.

`CITATION.cff` lists both manuscripts as unpublished. No acceptance, venue,
or DOI is claimed. The software and generated data use the MIT license;
this does not change rights in either manuscript or the LNCS class.
