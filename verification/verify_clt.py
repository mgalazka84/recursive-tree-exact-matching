"""Exact finite checks for the preferential cumulative-cost CLT proof.

Checks the nested-pair variance identity for size-centered history tolls,
the mean-preserving approximation, and the four-vertex resampling argument
for strictly positive variance. Finite checks do not prove the CLT.
Standard library only; no Monte Carlo tests of normality are used.
"""
from fractions import Fraction as F
from itertools import product
from math import prod
from pathlib import Path
import argparse
import json

from verify_preferential import history_weight, history_statistics, fringe_statistics, states, rational


def records(n):
    for prefix in product(*(range(v) for v in range(1, n))):
        parents = (-1,) + prefix
        children = [[] for _ in parents]
        for v in range(1, n):
            children[parents[v]].append(v)
        history = history_statistics(parents, children)
        sizes = fringe_statistics(children)[0]
        yield parents, children, history, sizes


def replay_root_signal(parents, full_trace=False):
    children, signal, total, old = [[]], [0], 0, [0]
    trace = [(0,)]
    for v in range(1, len(parents)):
        children.append([])
        children[parents[v]].append(v)
        new = states(children)
        changes = sum(a != b for a, b in zip(old, new))
        total += changes + 1 - int(old[0] != new[0])
        signal.append(new[0])
        trace.append(tuple(new))
        old = new
    return (tuple(signal), total, trace) if full_trace else (tuple(signal), total)


def check_gadget_replacements(parents, children, total):
    checked = 0
    for root in range(len(parents)):
        stack, vertices = [root], []
        while stack:
            v = stack.pop()
            vertices.append(v)
            stack.extend(children[v])
        vertices.sort()
        if len(vertices) != 4:
            continue
        r, a, b, c = vertices
        if parents[a] != r or parents[b] != r or parents[c] not in (r, a, b):
            continue
        original = int(parents[c] != r)
        root_signal, _, original_trace = replay_root_signal(parents, full_trace=True)
        interior = set(vertices)
        for new_parent in (r, a, b):
            altered = list(parents)
            altered[c] = new_parent
            signal, cost, altered_trace = replay_root_signal(altered, full_trace=True)
            assert signal == root_signal
            for before, after in zip(original_trace, altered_trace):
                assert all(before[v] == after[v] for v in range(len(before)) if v not in interior)
            assert cost - total == int(new_parent != r) - original
            checked += 1
    return checked


def verify(max_n=8):
    betas = tuple(map(F, ('0', '1/10', '1/4', '1/2', '1', '2', '10')))
    means = {b: {} for b in betas}
    etas = {b: {} for b in betas}
    variances, gadgets = [], []
    count = replacement_count = identity_count = 0
    for n in range(1, max_n + 1):
        trees = list(records(n))
        count += len(trees)
        for parents, children, history, nk in trees:
            replacement_count += check_gadget_replacements(parents, children, history[1])
        for beta in betas:
            weights = [history_weight(p, beta) for p, _, _, _ in trees]
            normalization = sum(weights)
            assert normalization == prod((beta.denominator + beta.numerator) * j - beta.numerator for j in range(1, n))
            means[beta][n] = F(sum(w * rec[2][0] for w, rec in zip(weights, trees)), normalization)
            g = means[beta]
            zeta, discount = 1 / (1 + beta), beta / (1 + beta)
            if n > 4:
                a_values = [rec[2][1] - n + 1 + rec[2][0] for rec in trees]
                ea = F(sum(w * a for w, a in zip(weights, a_values)), normalization)
                va = F(sum(w * a * a for w, a in zip(weights, a_values)), normalization) - ea ** 2
                eta_gadget = (1 + beta) / (2 + beta)
                v_gadget = 2 * (1 + 2 * beta) / (3 + 2 * beta) ** 2
                expected_four_fringes = (n - discount) * zeta / ((4 - discount) * (5 - discount))
                assert va >= eta_gadget * v_gadget * expected_four_fringes
            for cutoff in range(n):
                ef = variance = eta = F(0)
                for w, (_, _, history, nk) in zip(weights, trees):
                    root_score, cumulative, _, scores, _ = history
                    residual = sum((scores[k] - g[k] * nk[k] for k in range(cutoff + 1, n + 1)), F(0))
                    root_residual = root_score - g[n]
                    mass = F(w, normalization)
                    ef += mass * residual
                    variance += mass * residual ** 2
                    eta += mass * (2 * root_residual * residual - root_residual ** 2)
                    additive = cumulative - n + 1 + root_score
                    approximation = sum((scores[k] if k <= cutoff else g[k] * nk[k] for k in range(1, n + 1)), F(0))
                    assert residual == additive - approximation
                assert ef == 0, (beta, n, cutoff, 'residual mean')
                etas[beta][cutoff, n] = eta
                predicted = eta + (n - discount) * sum((
                    zeta * etas[beta].get((cutoff, k), F(0)) / ((k - discount) * (k + 1 - discount))
                    for k in range(1, n)), F(0))
                assert variance == predicted, (beta, n, cutoff, 'nested-pair identity')
                identity_count += 1
                if n == max_n:
                    variances.append({'beta': str(beta), 'n': n, 'cutoff': cutoff,
                                      'Var_centered_residual': rational(variance), 'root_eta': rational(eta)})
            if n == 4:
                group = [(w, history[1] - 3 + history[0]) for w, (parents, _, history, _) in zip(weights, trees)
                         if replay_root_signal(parents)[0] == (0, 1, 1, 1)]
                assert len(group) == 3
                group_weight = sum(w for w, _ in group)
                probability = F(group_weight, normalization)
                mean = F(sum(w * a for w, a in group), group_weight)
                variance = F(sum(w * a * a for w, a in group), group_weight) - mean ** 2
                assert probability == (1 + beta) / (2 + beta)
                assert variance == 2 * (1 + 2 * beta) / (3 + 2 * beta) ** 2
                lower = zeta * probability * variance / ((4 - discount) * (5 - discount))
                assert lower > 0
                gadgets.append({'beta': str(beta), 'signal_probability': rational(probability),
                                'conditional_variance': rational(variance), 'asymptotic_variance_lower_bound': rational(lower)})
        print(f'n={n}: centered variance identities and gadget replacements passed', flush=True)
    return {'max_n': max_n, 'beta_values': list(map(str, betas)), 'arithmetic': 'exact fractions',
            'distinct_labelled_trees': count, 'centered_variance_identities': identity_count,
            'gadget_replacements': replacement_count,
            'finite_residual_variances_at_max_n': variances, 'positive_variance_gadget': gadgets,
            'scope': 'Finite identities only; asymptotic CLT and variance bounds are proved in the manuscript.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-n', type=int, default=8)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = verify(args.max_n)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + '\n')
    print('All CLT-proof finite checks passed.')
