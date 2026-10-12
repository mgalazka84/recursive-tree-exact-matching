Repository note: verification/ paths below are relative to the repository
root. Compile the manuscript from this directory, or use make paper /
make paper-affine from the repository root. The downloadable Overleaf
archives also bundle their own verification scripts.

Exact Matching in Growing Recursive Trees with Constant Expected Update Cost
Marek Gałązka, Adam Mickiewicz University, Poznań, Poland

This is the original uniform-attachment article, restored without textual
changes from repository commit fd72ed52a1eea379ccde2fedcd7cfa2c409d2b00.
It has 12 LNCS pages in total, including the bibliography and all proofs.
The affine extension is a separate manuscript, affine_matching.tex.

Main document: recursive_matching.tex
Upload this ZIP to Overleaf and select that main document, or run:
  latexmk -pdf -interaction=nonstopmode -halt-on-error recursive_matching.tex
Alternatively run pdfLaTeX twice. Bibliography entries are embedded in the
source. llncs.cls is Springer LNCS v2.26; its original notices are retained.
The diagram is exact TikZ. No external raster figures are needed.

Contents of this paper:
- Theorem 1: exact expected update counts and constant expected recourse.
- Lemmas 2--3: correctness, matching representation, and implementation.
- Lemmas 4--6: root sensitivity, fringe means, and propagation generating functions.
- Theorem 7: limiting update distribution and moment convergence.
- Lemma 8: a fixed-size fringe variance bound with source attribution.
- Theorem 9: the quantitative L1 law for cumulative inspections.
- Theorem 10: the first-order maximum in probability.
- Figure 1 and Table 1: the local update example and uniform simulations.
No affine-attachment results or cumulative-cost CLT have been inserted here.

Verification and experiments (Python standard library, C++17):
  python3 verification/verify_extensions.py
  g++ -std=c++17 -O2 -Wall -Wextra -Wpedantic verification/matching_experiments.cpp -o verification/matching_experiments
  ./verification/matching_experiments verification
  python3 verification/summarize_experiments.py
Keep assertions enabled. The source includes a fallback copy of the table,
so it still compiles if verification/experiment_table.tex was omitted.

The uniform data contain 32 independent trajectories with seeds 0--31 and
final size 1,000,000. Checkpoints on the same trajectory are dependent.
The finite checks and simulations supplement, rather than replace, proofs.

Public implementation and reproduction entry point:
  https://github.com/mgalazka84/recursive-tree-exact-matching
  python3 reproduce.py
The original uniform-experiment snapshot cited in the article is
0e4c10c707699f61ca2e83cd28de8ca4683ef03d.
The repository retains common algorithm code and separates the manuscripts
into paper/ (uniform) and paper_affine/ (affine preferential attachment).
