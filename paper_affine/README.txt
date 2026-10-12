Repository note: verification/ paths below are relative to the repository
root. Compile the manuscript from this directory, or use make paper /
make paper-affine from the repository root. The downloadable Overleaf
archives also bundle their own verification scripts.

Exact Matching under Affine Preferential Attachment:
Update Laws, Gaussian Fluctuations, and Extremes
Marek Gałązka, Adam Mickiewicz University, Poznań, Poland

Standalone affine-attachment article: 10 main-text pages, 1 bibliography
page, and 7 appendix pages (18 pages total in standard LNCS format).
Main document: affine_matching.tex

Upload this ZIP to Overleaf and select affine_matching.tex, or run:
  latexmk -pdf -interaction=nonstopmode -halt-on-error affine_matching.tex
Alternatively run pdfLaTeX twice. Bibliography entries are embedded.
llncs.cls is Springer LNCS v2.26 with its original notices retained.
Figure 1 is an exact TikZ rendering of the three four-vertex histories used
to prove positive variance. It has no external image dependency.

Model: parent weight 1 + beta * number_of_children, fixed beta > 0.
The separate seed-edge model uses degree + delta, delta > -1, and maps to
beta = 1/(1+delta) in its two components.

Results and proof locations:
- Lemma 1: recalled deterministic matching algorithm and cost conventions.
- Theorem 2: constant expected update/recourse and factorial tail bound;
  proof in Appendix A.
- Theorem 3: limiting update distribution, moment convergence, sharp tail;
  proof and certified rational truncations in Appendix B.
- Theorem 4: almost-sure bounded window for the maximum, with Lambert-W
  centering and second-order expansion; proof in Appendix C.
- Theorem 5: linear cumulative mean and explicit constant correction;
  proof in Appendix D.
- Theorem 6: full cumulative-cost CLT, linear positive variance, and the
  exact asymptotic absolute mean error; complete proof in Section 4.
- Corollary 7: degree-based attachment, including E U_n <= 3/2 for ordinary
  degree attachment; component comparisons in Appendix E and Section 4.
- Lemma 8: fringe means, a cited uniform size-count variance estimate,
  and conditional independence of disjoint interiors.
- Lemma 9: variance bound for size-centered additive tolls.

Self-contained scope:
All definitions, statements and proofs needed for the affine results are
included. The classical matching recurrence is attributed to the existing
literature, and the implementation lemma is proved directly. The text does
not cite or refer to the separate uniform-attachment manuscript. The beta=0
CLT is not claimed. No closed-form variance constant, Gumbel limit, or
beta-uniform error bound is claimed.

Exact finite verification (Python 3, standard library only):
  python3 verification/verify_preferential.py
  python3 verification/verify_clt.py
The second script imports helpers from the first; keep them together.
Saved reference data are verification/preferential_exact_checks.json and
verification/clt_exact_checks.json. The scripts are identical to the shared
repository versions. They enumerate 5,914 trees and 46,233 insertions for
seven beta values and four degree offsets. The CLT checker additionally
checks 252 centered variance identities and 3,474 replacements preserving
every exterior state at every arrival prefix. Finite checks do not prove
asymptotic normality; the mathematical proof is included in the manuscript.

Public implementation and one-command reproduction:
  https://github.com/mgalazka84/recursive-tree-exact-matching
  python3 reproduce.py --checks-only
The repository has paper/ for the uniform article and paper_affine/ for
this article, sharing include/, verification/, and data/.
