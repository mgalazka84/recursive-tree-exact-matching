# Exact matching in growing random trees

Reproducibility code for **Exact Matching in Growing Random Trees with Constant Expected Update Cost**, by **Marek Gałązka**, Adam Mickiewicz University, Poznań.

Public repository: https://github.com/mgalazka84/recursive-tree-exact-matching.
The manuscript retains snapshot `0e4c10c707699f61ca2e83cd28de8ca4683ef03d` for the original uniform-attachment experiments. The current repository also includes the affine preferential-attachment extension and its exact checker.

The code maintains an actual maximum-cardinality matching while vertices arrive as leaves. It includes a reusable C++ implementation, constant-time size/mate/edge queries, independent finite checks, and all data and scripts behind the manuscript's simulation table. Mathematical proofs are in the manuscript; computational checks supplement them.

## Quick start

Requirements: **Python 3.8 or later**, a **C++17 compiler** (such as `g++`), and no external Python packages. LaTeX is needed only to build the paper.

The supplied reproduction was checked with Python 3.12.14 and GNU C++ 13.3.0.

From the repository root:

```sh
python3 reproduce.py
```

This compiles the implementation, checks the query API against brute-force edge-subset optima, runs the Python finite checkers (including preferential attachment), enumerates all 46,234 recursive trees of sizes 1--9 using the complete matching implementation, and regenerates 32 growth trajectories to one million vertices. New raw data and aggregate results are written to `results/`, compared with the saved data in `data/`, and the manuscript table is refreshed.

For the finite checks alone:

```sh
python3 reproduce.py --checks-only
```

The finite Python checker also retains its original 12 smaller supporting runs to 100,000 vertices. These are separate from the 32 runs used in the article's table.

Optional Make targets:

```sh
make demo       # the exact example in Figure 1
make checks     # query/brute-force and mathematical finite checks
make reproduce  # all experiments and reference-data comparisons
make paper     # compile the LNCS source; requires latexmk and TikZ
```

Keep C++ assertions enabled; do not compile verification executables with `-DNDEBUG`. If another compiler is preferred, set `CXX=clang++` before running the reproduction script.

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

Correctness holds for **every** sequence of leaf insertions. The paper proves detailed exact and limit laws for **independent uniform parent choices**, and constant expected maintenance under **affine preferential attachment**. Expectations average over the entire growth history, rather than conditioning on every possible realized tree. Edge weights and deletions are outside the guarantees. The parent is supplied to the algorithm; generating a random parent is not included in its update cost.

For the one-root model with parent weight `1 + beta * number_of_children`, `beta >= 0`, let `Z_n = (1+beta)n - beta`. Theorem 11 proves `E K_n <= n/Z_n`, `E U_n <= 1+(n-1)/Z_n`, and expected edge recourse at most `2+n/Z_n <= 3`. For each fixed positive beta it also gives the uniform tail `P(K_n >= k) <= exp(-k log k + O_beta(k))` and the corresponding maximum upper bound.

For the distinct seed-edge model with weights `degree + delta`, `delta > -1`, Corollary 12 gives `E U_n <= 1+(1+delta)(n-1)/((2+delta)n-2)`. In particular, ordinary degree-proportional attachment has `E U_n <= 3/2`. The proof excludes the overall root when counting changed states, so it does not conflate the two root conventions. The uniform limiting law, total-cost concentration rate, and sharp maximum lower bound have not been proved for these preferential models.

The reference code uses `std::vector` storage; an optional capacity reservation is used for fixed-horizon simulations. Without reservation, allocator effects have amortized rather than per-insertion guarantees. The paper's experimental counters measure state propagation and recourse, not allocator cost or wall-clock time. `verify_invariants()` is an independent O(n) diagnostic and is excluded from propagation counts.

## Files and mathematical results

| File | Role | Manuscript connection |
| --- | --- | --- |
| `include/recursive_tree_matching.hpp` | Complete algorithm and public queries | Lemmas 2--3; algorithmic part of Theorem 1 |
| `examples/demo.cpp` | Concrete insertion and its maximum matchings | Figure 1 |
| `tests/test_matching.cpp` | Edge-subset optimum enumeration and every query on all 874 trees of sizes 1--7 | Independent correctness checks |
| `verification/matching_experiments.cpp` | Full matching checks on 46,234 trees; exact sums; 32 uniform growth runs | Theorem 1; experimental checks related to Theorems 7, 9, 10 |
| `verification/verify_extensions.py` | Independently replayed subtree histories, exact covariance and fringe-path checks, finite PGFs | Lemmas 4, 6, 8; history identity in Theorem 9; path lower bound in Theorem 10 |
| `verification/verify_preferential.py` | Exact weighted enumeration for seven beta values and four degree-offset values; formal series through degree 24 | Theorem 11, Corollary 12, Appendix D |
| `verification/summarize_experiments.py` | Aggregate raw output and regenerate the simulation table | Section 10, Table 1 |
| `reproduce.py` | Compile, run checks, regenerate outputs and compare saved reference data | Complete reproduction entry point |
| `paper/recursive_matching.tex` | LNCS manuscript: 12 main-text pages, 1 references page, and 4 appendix pages; all proofs and the exact TikZ diagram | Definitions, statements, and proofs |

The root-change probability checked in exact enumeration is `1/n`, equivalent to the expected root sensitivity being one. The finite-PGF checker compares enumerated update distributions with independently expanded generating functions. These checks are not formal verification or asymptotic proofs.

## Exact preferential-attachment checks

```sh
python3 verification/verify_preferential.py --output results/preferential_exact_checks.json
```

The checker enumerates all **5,914** increasing labelled trees with old size 1--8 and all **46,233** possible insertions. It recomputes the full state vector after each insertion, independently of the sensitivity recurrence. Integer history weights and exact rational arithmetic cover beta values `0, 1/10, 1/4, 1/2, 1, 2, 10` and ordinary-degree offsets `-1/2, 0, 1, 3`.

It checks normalization, weighted root sensitivities, the finite affine fringe identity, three marked EGFs (`s=1,2,3`), and the mean bounds. A separate formal-series computation verifies coefficient domination through degree 24. Saved exact results are in `data/preferential_exact_checks.json`; the one-command reproduction compares them exactly. These are finite checks, not proofs of asymptotic tail estimates. No preferential-attachment Monte Carlo results are claimed.

## Saved data and simulation conventions

The initial tree is one root. At old size `n`, the new vertex chooses a parent uniformly from `0,...,n-1`. The deterministic RNG is the standard `std::mt19937_64`, with seeds **0--31**. Explicit rejection sampling avoids implementation-specific `uniform_int_distribution` sequences.

Each of the same 32 trajectories is observed at sizes `10^3`, `10^4`, `10^5`, and `10^6`. Updates within a trajectory, and observations at different sizes of that trajectory, are dependent. Table 1's SD is the **sample standard deviation across the 32 runs**. It is not an update-level standard error. No independence-based confidence interval is attached to pooled update moments.

| File in `data/` | Contents |
| --- | --- |
| `matching_exact_checks.csv` | Integer exact sums for old sizes n=1,...,8 |
| `matching_runs.csv` | One row per seed/checkpoint: S/N, W, W/ell, mean recourse, and late-update moments |
| `matching_late_histogram.csv` | Pooled K histogram from the last 500,000 updates in each run (16 million observations) |
| `experiment_summary.json` | Aggregates, seeds, RNG conventions, and limiting constants |
| `experiment_table.tex` | Table generated directly from the raw CSV files |

`verification/verification_results.json` records the separate original Python finite checks and supporting simulations. To regenerate the aggregate table without rerunning the large simulations, use:

```sh
python3 verification/summarize_experiments.py data
```

At N=10^6, the mean S/N is **1.9998017**, mean recourse is **1.307612**, and mean W/ell is **2.126309** (W ranges from 10 to 13). Pooled late-update variance is **1.192112**, compared with the limiting value **1.192695**.

The observed recourse does not improve the proved upper bound three. The largest-update theorem is a first-order limit with leading constant one; its asymptotic constant is still a poor numerical approximation at these sizes. The data establish neither a convergence rate nor monotonicity. The total-cost rate is an upper bound without an optimality claim; a central limit theorem is not established in this manuscript.

## Building and citing the manuscript

```sh
cd paper
latexmk -pdf -interaction=nonstopmode -halt-on-error recursive_matching.tex
```

The bibliography is embedded in the source. `paper/llncs.cls` retains Springer's original notices and is distributed under CC BY 4.0; see `THIRD_PARTY_NOTICES.md`. The tree figure uses TikZ, so no external raster image is required. The generated table is included from `paper/verification/experiment_table.tex` when available. The main source also embeds an identical fallback copy of Table 1, allowing compilation when that extra file was omitted from an Overleaf upload.

`CITATION.cff` gives the author and manuscript title. It does not claim an acceptance, publication venue, or DOI. The original software and generated data use the MIT license; that license does not change the rights in the manuscript or third-party LNCS class.
