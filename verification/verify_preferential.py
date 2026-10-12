"""Exact checks of the affine-attachment extension (standard library only).

Enumerate increasing labelled trees and recompute states after every possible
insertion, independently of the sensitivity recurrence. Compare weighted sums
with the formal differential equations and the affine fringe identity.
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
    marked = {}
    for s in (1, 2, 3):
        exp = phi(h, F(0), F(s), degree)
        exp[0] -= 1
        marked[s] = [x / s for x in exp]
    assert all(x <= y for x, y in zip(marked[1], tree))
    assert beta or marked[1] == tree
    return tree, marked


def states(children):
    answer = [0] * len(children)
    for v in reversed(range(len(children))):
        answer[v] = int(any(answer[u] == 0 for u in children[v]))
    return answer


def tree_records(n):
    for prefix in product(*(range(i) for i in range(1, n))):
        parents = (-1,) + prefix
        children = [[] for _ in range(n)]
        depth = [0] * n
        for v in range(1, n):
            children[parents[v]].append(v)
            depth[v] = depth[parents[v]] + 1
        old = states(children)
        costs = []
        for p in range(n):
            augmented = [list(xs) for xs in children] + [[]]
            augmented[p].append(n)
            new = states(augmented)
            changed = sum(x != y for x, y in zip(old, new))
            root = int(old[0] != new[0])
            costs.append((changed, changed + 1 - root, root))
        yield parents, [len(xs) for xs in children], depth, costs


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
    rows, degree_rows = [], []
    tree_count = update_count = 0
    for n in range(1, max_n + 1):
        sums = {b: [0] * 8 for b in betas}
        degree_sums = {d: [0] * 3 for d in deltas}
        for parents, children, depths, costs in tree_records(n):
            tree_count += 1
            update_count += n
            for b in betas:
                a, den = b.numerator, b.denominator
                w = history_weight(parents, b)
                total = sums[b]
                total[0] += w
                for p, (k, u, root) in enumerate(costs):
                    mass = w * (den + a * children[p])
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
            rows.append({"beta": str(b), "n": n, "E_B": rational(eb), "E_K": rational(ek), "E_U": rational(eu)})
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
    return {"max_n": max_n, "series_degree": series_degree,
            "distinct_labelled_trees": tree_count, "distinct_insertions": update_count,
            "beta_values": list(map(str, betas)), "delta_values": list(map(str, deltas)),
            "arithmetic": "exact fractions and integer history weights",
            "checks": ["normalization", "independent state recomputation", "weighted root sensitivity",
                       "three marked generating functions", "affine fringe identity", "per-update bounds",
                       "formal coefficient domination", "ordinary-degree seed-edge corollary"],
            "affine_outdegree": rows, "ordinary_degree": degree_rows}


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
