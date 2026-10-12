"""Exact checks of the affine-attachment extension (standard library only).

Enumerate increasing labelled trees and recompute states after every possible
insertion, independently of the sensitivity recurrence. Compare weighted sums
with the formal differential equations and the affine fringe identity.
Also check every finite update-tail probability, fringe-path thinning,
and the pathwise comparison with the two components of a seed-edge tree.
Rational truncations give certified intervals for the limiting law.
Replay all relative-label subtree histories to check cumulative means,
fixed-size history-score variances, and last-insertion witnesses for maxima.
These finite checks supplement, and do not replace, the proofs.
"""
from fractions import Fraction as F
from itertools import product
from math import factorial, prod
from pathlib import Path
import argparse
import json


def mul(a, b, degree):
    out = [F(0)] * (degree + 1)
    for i, x in enumerate(a[:degree + 1]):
        if x:
            for j, y in enumerate(b[:degree + 1 - i]):
                out[i + j] += x * y
    return out


def phi(a, beta, gamma, degree):
    """(1-beta*a)^(-gamma/beta), or exp(gamma*a) if beta=0."""
    assert not a[0]
    out = [F(0)] * (degree + 1)
    power = [F(1)] + [F(0)] * degree
    coefficient = F(1)
    for j in range(degree + 1):
        for i, x in enumerate(power):
            out[i] += coefficient * x
        power = mul(power, a, degree)
        coefficient *= (gamma + beta * j) / (j + 1)
    return out


def series(beta, degree):
    tree = [F(0)] * (degree + 1)
    weight = F(1)
    for n in range(1, degree + 1):
        tree[n] = weight / factorial(n)
        weight *= (1 + beta) * n - beta
    zero = [F(0)] * (degree + 1)
    for j in range(degree):
        one = [tree[i] - zero[i] for i in range(j + 1)]
        zero[j + 1] = phi(one, beta, F(1), j)[j] / (j + 1)
    one = [x - y for x, y in zip(tree, zero)]
    hprime = phi(one, beta, 1 + beta, degree - 1)
    h = [F(0)] + [x / (j + 1) for j, x in enumerate(hprime)]
    marked = {0: h}
    for s in (1, 2, 3):
        exp = phi(h, F(0), F(s), degree)
        exp[0] -= 1
        marked[s] = [x / s for x in exp]
    assert all(x <= y for x, y in zip(marked[1], tree))
    assert beta or marked[1] == tree
    return tree, marked


def fringe_statistics(children):
    n = len(children)
    sizes, paths = [1] * n, [True] * n
    counts, path_counts = [0] * (n + 1), [0] * (n + 1)
    for v in reversed(range(n)):
        sizes[v] = 1 + sum(sizes[u] for u in children[v])
        paths[v] = len(children[v]) <= 1 and all(paths[u] for u in children[v])
        counts[sizes[v]] += 1
        path_counts[sizes[v]] += int(paths[v])
    return counts, path_counts


def fringe_average(beta, n, value):
    q = beta / (1 + beta)
    total = sum((n - q) * (1 - q) / ((j - q) * (j + 1 - q)) * value(j)
                for j in range(1, n)) + value(n)
    return total / ((1 + beta) * n - beta)


def tail_series(h, max_k, degree):
    power = [F(1)] + [F(0)] * degree
    out = {}
    for k in range(1, max_k + 1):
        power = mul(power, h, degree)
        out[k] = [x / factorial(k) for x in power]
    return out


def certified_limit_interval(beta, coefficients, upper_cap):
    """Nonnegative truncation plus the proved coefficient-majorant remainder."""
    zeta = 1 / (1 + beta)
    degree = len(coefficients) - 1
    d, zpower, lower = F(1), F(1), F(0)
    for j in range(1, degree + 1):
        d *= F(j) / (j + zeta)
        zpower *= zeta
        lower += d * zpower * coefficients[j]
    error = zeta ** 2 / (degree + zeta)
    upper = min(upper_cap, lower + error)
    assert 0 <= lower <= upper
    return {"lower": rational(lower), "upper": rational(upper)}


def states(children):
    answer = [0] * len(children)
    for v in reversed(range(len(children))):
        answer[v] = int(any(answer[u] == 0 for u in children[v]))
    return answer


def history_statistics(parents, children):
    """Recompute every prefix of every fringe, without propagation code."""
    n = len(parents)
    scores = [0] * n
    by_size = [0] * (n + 1)
    witnesses = [[0] * n for _ in range(n + 1)]
    full_cost = full_maximum = root_witness = 0
    for root in range(n):
        vertices, stack = [], [root]
        while stack:
            v = stack.pop()
            vertices.append(v)
            stack.extend(children[v])
        vertices.sort()
        old = {root: 0}
        last_witness = 0
        for j in range(1, len(vertices)):
            present = set(vertices[:j + 1])
            new = {}
            for v in reversed(vertices[:j + 1]):
                new[v] = int(any(new[u] == 0 for u in children[v] if u in present))
            changed = sum(old[v] != new[v] for v in old)
            root_changed = int(old[root] != new[root])
            scores[root] += root_changed
            last_witness = changed if root_changed else 0
            if root == 0:
                cost = changed + 1 - root_changed
                full_cost += cost
                full_maximum = max(full_maximum, cost)
            old = new
        size = len(vertices)
        by_size[size] += scores[root]
        if last_witness:
            witnesses[size][last_witness] += 1
            assert full_maximum >= last_witness
        if root == 0:
            root_witness = last_witness
    assert full_cost == n - 1 + sum(scores) - scores[0]
    return scores[0], full_cost, root_witness, by_size, witnesses


def tree_records(n):
    for prefix in product(*(range(i) for i in range(1, n))):
        parents = (-1,) + prefix
        children = [[] for _ in range(n)]
        depth = [0] * n
        for v in range(1, n):
            children[parents[v]].append(v)
            depth[v] = depth[parents[v]] + 1
        old = states(children)
        if n >= 2:
            forest = [list(xs) for xs in children]
            forest[0].remove(1)
            forest_old = states(forest)
            component = [0, 1] + [0] * (n - 2)
            for v in range(2, n):
                component[v] = component[parents[v]]
        costs = []
        for p in range(n):
            augmented = [list(xs) for xs in children] + [[]]
            augmented[p].append(n)
            new = states(augmented)
            changed = sum(x != y for x, y in zip(old, new))
            root = int(old[0] != new[0])
            costs.append((changed, changed + 1 - root, root))
            if n >= 2:
                forest_augmented = [list(xs) for xs in forest] + [[]]
                forest_augmented[p].append(n)
                forest_new = states(forest_augmented)
                kf = sum(a != b for a, b in zip(forest_old, forest_new))
                rf = int(forest_old[component[p]] != forest_new[component[p]])
                uf = kf + 1 - rf
                assert abs(changed - kf) <= rf
                assert changed + 1 - root == uf + (rf if component[p] == 1 else 0)
        yield (parents, [len(xs) for xs in children], depth, costs,
               fringe_statistics(children), history_statistics(parents, children))


def history_weight(parents, parameter, ordinary_degree=False):
    numerator, denominator = parameter.numerator, parameter.denominator
    count = [0] * len(parents)
    value = 1
    for v in range(1, len(parents)):
        p = parents[v]
        if not ordinary_degree:
            value *= denominator + numerator * count[p]
        elif v > 1:  # the first edge is the deterministic seed
            degree = count[p] + int(p != 0)
            value *= denominator * degree + numerator
        count[p] += 1
    return value


def rational(x):
    return {"exact": str(x), "decimal": float(x)}


def verify(max_n=8, series_degree=24):
    betas = tuple(map(F, ("0", "1/10", "1/4", "1/2", "1", "2", "10")))
    deltas = tuple(map(F, ("-1/2", "0", "1", "3")))
    formal = {b: series(b, series_degree) for b in betas}
    tails = {b: tail_series(formal[b][1][0], max_n, series_degree) for b in betas}
    rows, degree_rows, fringe_rows, witness_rows = [], [], [], []
    root_moments = {b: {} for b in betas}
    tree_count = update_count = 0
    for n in range(1, max_n + 1):
        sums = {b: [0] * 8 for b in betas}
        degree_sums = {d: [0] * 3 for d in deltas}
        path_sums = {b: [[0] * 4 for _ in range(n + 1)] for b in betas}
        histograms = {b: [0] * (n + 1) for b in betas}
        histories = {b: [0] * 4 for b in betas}
        last_histograms = {b: [0] * n for b in betas}
        score_sums = {b: [[0, 0] for _ in range(n + 1)] for b in betas}
        witness_sums = {b: [[[0, 0] for _ in range(n)] for _ in range(n + 1)] for b in betas}
        for parents, children, depths, costs, (nk, xk), history in tree_records(n):
            tree_count += 1
            update_count += n
            root_score, cumulative, last_witness, yk, witnesses = history
            for b in betas:
                a, den = b.numerator, b.denominator
                w = history_weight(parents, b)
                total = sums[b]
                total[0] += w
                for j, x in enumerate((root_score, root_score ** 2, cumulative, cumulative ** 2)):
                    histories[b][j] += w * x
                last_histograms[b][last_witness] += w
                for k in range(1, n + 1):
                    for j, x in enumerate((nk[k], nk[k] ** 2, xk[k], xk[k] ** 2)):
                        path_sums[b][k][j] += w * x
                    score_sums[b][k][0] += w * yk[k]
                    score_sums[b][k][1] += w * yk[k] ** 2
                    if n == max_n:
                        for r in range(1, k):
                            x = witnesses[k][r]
                            witness_sums[b][k][r][0] += w * x
                            witness_sums[b][k][r][1] += w * x ** 2
                for p, (k, u, root) in enumerate(costs):
                    mass = w * (den + a * children[p])
                    histograms[b][k] += mass
                    total[1] += mass * k
                    total[2] += mass * u
                    for offset, s in enumerate((1, 2, 3)):
                        total[3 + offset] += mass * root * s ** depths[p]
                    total[6] += mass * (2 ** k - 1)
                    total[7] += mass * ((3 ** k - 1) // 2)
            if n >= 2:
                for d in deltas:
                    w = history_weight(parents, d, True)
                    degree_sums[d][0] += w
                    for p, (k, u, _) in enumerate(costs):
                        mass = w * (d.denominator * (children[p] + int(p != 0)) + d.numerator)
                        degree_sums[d][1] += mass * k
                        degree_sums[d][2] += mass * u
        for b, total in sums.items():
            a, den = b.numerator, b.denominator
            normalization = prod((den + a) * j - a for j in range(1, n))
            zn_integer = (den + a) * n - a
            zn = F(zn_integer, den)
            assert total[0] == normalization
            ek, eu = (F(total[i], normalization * zn_integer) for i in (1, 2))
            assert ek <= F(n) / zn
            assert eu <= 1 + F(n - 1) / zn
            tree, marked = formal[b]
            eq, eq2, es, es2 = [F(x, normalization) for x in histories[b]]
            root_moments[b][n] = (eq, eq2)
            zeta = 1 / (1 + b)
            d, zpower = F(1), F(1)
            predicted_eu, predicted_es, predicted_eq = F(1), F(n - 1), F(0)
            for j in range(1, n):
                d *= F(j) / (j + zeta)
                zpower *= zeta
                a_j = marked[1][j] * zpower
                predicted_eu += d * a_j
                predicted_es += (n - 1 - j) * d * a_j
                predicted_eq += marked[1][j] / tree[j] / ((1 + b) * j - b)
            assert eu == predicted_eu, (b, n, "exact inspection prefix")
            assert es == predicted_es, (b, n, "exact cumulative mean")
            assert eq == predicted_eq, (b, n, "root-history mean")
            for r in range(1, n):
                probability = F(last_histograms[b][r], normalization)
                assert probability == tails[b][r][n - 1] / (n * tree[n]), (b, n, r, "integrated witness EGF")
            eb = F(total[3], normalization * den)
            assert eb <= 1
            if b == 0:
                assert eb == ek == 1 and eu == 2 - F(1, n)
            for offset, s in enumerate((1, 2, 3)):
                ebs = F(total[3 + offset], normalization * den)
                assert ebs == marked[s][n] / tree[n], (b, n, s, "marked EGF")
                q = b / (1 + b)
                fringe_sum = sum(
                    (n - q) * (1 - q) / ((j - q) * (j + 1 - q)) * marked[s][j] / tree[j]
                    for j in range(1, n)
                ) + ebs
                exact_f = ek if s == 1 else F(total[5 + offset], normalization * zn_integer)
                assert exact_f == fringe_sum / zn, (b, n, s, "fringe identity")
            for k in range(1, n + 1):
                enumerated_tail = F(sum(histograms[b][k:]), normalization * zn_integer)
                predicted_tail = fringe_average(b, n, lambda j: tails[b][k][j] / tree[j])
                assert enumerated_tail == predicted_tail, (b, n, k, "finite update tail")
                en, en2, ex, ex2 = [F(x, normalization) for x in path_sums[b][k]]
                pi = F(1) / prod((1 + b) * i - b for i in range(1, k))
                vn, vx = en2 - en ** 2, ex2 - ex ** 2
                assert ex == pi * en
                assert vx == pi ** 2 * vn + pi * (1 - pi) * en
                ey, ey2 = [F(x, normalization) for x in score_sums[b][k]]
                g, g2 = root_moments[b][k]
                assert ey == en * g
                assert ey2 - ey ** 2 == en * (g2 - g ** 2) + vn * g ** 2
                if k < n:
                    q = b / (1 + b)
                    assert en == (n - q) * (1 - q) / ((k - q) * (k + 1 - q))
                if n == max_n:
                    fringe_rows.append({"beta": str(b), "n": n, "k": k,
                                        "path_probability": rational(pi),
                                        "E_paths": rational(ex), "Var_paths": rational(vx)})
                    for r in range(1, k):
                        probability = tails[b][r][k - 1] / (k * tree[k])
                        ew, ew2 = [F(x, normalization) for x in witness_sums[b][k][r]]
                        vw = ew2 - ew ** 2
                        assert ew == en * probability
                        assert vw == probability ** 2 * vn + probability * (1 - probability) * en
                        witness_rows.append({"beta": str(b), "n": n, "fringe_size": k,
                                             "root_sensitive_flips": r, "probability": rational(probability),
                                             "E_witnesses": rational(ew), "Var_witnesses": rational(vw)})
            rows.append({"beta": str(b), "n": n, "E_B": rational(eb), "E_K": rational(ek), "E_U": rational(eu),
                         "E_root_history": rational(eq), "E_root_history_squared": rational(eq2),
                         "E_S": rational(es), "Var_S": rational(es2 - es ** 2)})
        for d, total in degree_sums.items():
            if n < 2:
                continue
            den, a = d.denominator, d.numerator
            normalization = prod((2 * den + a) * j - 2 * den for j in range(2, n))
            zn = (2 + d) * n - 2
            assert total[0] == normalization
            ek, eu = (F(total[i], normalization * ((2 * den + a) * n - 2 * den)) for i in (1, 2))
            bound = 1 + (n - 1) * (1 + d) / zn
            assert eu <= bound and ek <= eu
            if d == 0:
                assert eu <= F(3, 2)
            degree_rows.append({"delta": str(d), "n": n, "E_K": rational(ek), "E_U": rational(eu), "U_bound": rational(bound)})
        print(f"n={n}: all affine and ordinary-degree checks passed", flush=True)
    intervals = []
    for b in betas:
        zeta = 1 / (1 + b)
        intervals.append({"beta": str(b),
                          "mean_K": certified_limit_interval(b, formal[b][1][1], zeta),
                          "tail_probabilities": {
                              str(k): certified_limit_interval(b, tails[b][k], F(1))
                              for k in range(1, max_n + 1)}})
    return {"max_n": max_n, "series_degree": series_degree,
            "distinct_labelled_trees": tree_count, "distinct_insertions": update_count,
            "beta_values": list(map(str, betas)), "delta_values": list(map(str, deltas)),
            "arithmetic": "exact fractions and integer history weights",
            "checks": ["normalization", "independent state recomputation", "weighted root sensitivity",
                       "three marked generating functions", "affine fringe identity", "per-update bounds",
                       "formal coefficient domination", "ordinary-degree seed-edge corollary",
                       "all finite update-tail coefficients", "fringe-path mean and variance thinning",
                       "pointwise seed-edge forest comparison", "certified limiting-law intervals",
                       "relative-label history identity", "exact cumulative mean and inspection prefix",
                       "compound fringe-history variance", "integrated last-insertion witness EGF",
                       "last-insertion witness thinning and pathwise maximum lower bound"],
            "affine_outdegree": rows, "ordinary_degree": degree_rows,
            "fringe_paths_at_max_n": fringe_rows,
            "last_insertion_witnesses_at_max_n": witness_rows,
            "limiting_law_intervals": intervals,
            "interval_remainder_bound": "zeta^2/(series_degree+zeta), zeta=1/(1+beta)"}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-n", type=int, default=8)
    parser.add_argument("--series-degree", type=int, default=24)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.series_degree < args.max_n:
        parser.error("series degree must be at least max-n")
    result = verify(args.max_n, args.series_degree)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"Verified {result['distinct_labelled_trees']:,} trees and {result['distinct_insertions']:,} insertions.")
