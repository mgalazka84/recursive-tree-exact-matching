"""Independent finite checks of the three new mathematical arguments."""
from collections import Counter
from fractions import Fraction as Q
from itertools import product
from math import factorial
from pathlib import Path
import json
import random


def harmonic(n):
    return sum((Q(1, j) for j in range(1, n+1)), Q(0))


def replay(parents):
    n = len(parents)
    children = [[] for _ in parents]
    states = [0]*n
    flips = [0]*n
    costs = []
    changes = []
    for x in range(1, n):
        children[parents[x]].append(x)
        v = parents[x]
        u = k = 0
        while True:
            u += 1
            new = int(any(states[c] == 0 for c in children[v]))
            if new == states[v]:
                break
            states[v] = new
            flips[v] += 1
            k += 1
            if v == 0:
                break
            v = parents[v]
        costs.append(u)
        changes.append(k)
    return children, flips, costs, changes


def independent_subtree_history(v, children, parents):
    stack = [v]
    vertices = []
    while stack:
        w = stack.pop()
        vertices.append(w)
        stack.extend(children[w])
    vertices.sort()
    old = count = 0
    for i in range(1, len(vertices)):
        present = set(vertices[:i+1])
        states = {}
        for w in reversed(vertices[:i+1]):
            states[w] = int(any(states[c] == 0 for c in children[w]
                                if c in present))
        count += states[v] != old
        old = states[v]
    path = all(len(children[w]) <= 1 for w in vertices)
    return len(vertices), count, path


records = []
root_moments = {}
last_histograms = {}
checked = 0
for n in range(1, 10):
    trees = factorial(n-1)
    hist = Counter()
    root_sum = root_sq = 0
    fixed = {k: [0, 0, 0, 0] for k in range(1, (n-1)//2+1)}
    for seq in product(*(range(i) for i in range(1, n))):
        parents = [-1] + list(seq)
        children, flips, costs, changes = replay(parents)
        summaries = [independent_subtree_history(v, children, parents)
                     for v in range(n)]
        q = [x[1] for x in summaries]
        assert q == flips
        assert sum(changes) == sum(q)
        assert sum(costs) == n-1 + sum(q) - q[0]
        longest_path = max(x[0] for x in summaries if x[2])
        assert max(costs, default=0) >= longest_path - 1
        root_sum += q[0]
        root_sq += q[0]**2
        if n > 1:
            hist[changes[-1]] += 1
        for k, stats in fixed.items():
            y = sum(x[1] for x in summaries if x[0] == k)
            p = sum(x[2] for x in summaries if x[0] == k)
            stats[0] += y
            stats[1] += y*y
            stats[2] += p
            stats[3] += p*p
        checked += 1
    mean = Q(root_sum, trees)
    second = Q(root_sq, trees)
    assert mean == harmonic(n-1)
    root_moments[n] = (mean, second)
    for k, (ys, yss, xs, xss) in fixed.items():
        beta = Q(3*k+2, k*(k+1)**2*(2*k+1))
        eq, eq2 = root_moments[k]
        ey = Q(ys, trees)
        vy = Q(yss, trees) - ey*ey
        assert ey == n*eq/Q(k*(k+1))
        assert vy == n*eq2/Q(k*(k+1)) - n*beta*eq*eq
        pp = Q(1, factorial(k-1))
        ex = Q(xs, trees)
        vx = Q(xss, trees)-ex*ex
        assert ex == Q(n, (k+1)*factorial(k))
        assert vx == ex - n*beta*pp*pp
    if n > 1:
        last_histograms[n-1] = {k: Q(v, trees) for k,v in hist.items()}
        assert sum(k*v for k,v in last_histograms[n-1].items()) == 1
        assert sum(2**k*v for k,v in last_histograms[n-1].items()) == 3-Q(1,n-1)
    records.append({'size':n, 'trees':trees,
                    'mean_root_flips':str(mean),
                    'second_root_moment':str(second)})
    print('exhaustive', n, trees, 'OK', flush=True)

# Obtain finite probability generating functions directly from the EGF,
# independently of the state-propagation enumeration above.
limit = 8
r = [Q(0)] + [Q(1,j) for j in range(1,limit+1)]
power = [Q(1)] + [Q(0)]*limit
b = [[Q(0)]*(limit+1) for _ in range(limit+1)]
falling = [Q(1)]
for j in range(1,limit+1):
    power = [sum((power[h]*r[k-h] for h in range(k+1)),Q(0))
             for k in range(limit+1)]
    if j > 1:
        new = [Q(0)]*(len(falling)+1)
        for h,c in enumerate(falling):
            new[h] += (1-j)*c
            new[h+1] += c
        falling = new
    for k in range(1,limit+1):
        for h,c in enumerate(falling):
            b[k][h] += power[k]*c/factorial(j)
for n, distribution in last_histograms.items():
    f = [sum((b[k][h]/(k+1) for k in range(1,n)),Q(0))+b[n][h]
         for h in range(limit+1)]
    generating = [Q(1)-f[0]]+[f[h-1]-f[h] for h in range(1,limit+1)]
    for k in range(n+1):
        expected = distribution.get(k,Q(0))
        assert generating[k] == expected
print('finite PGFs for n=1..8 agree exactly with exhaustive enumeration', flush=True)


def simulate(n, seed):
    rng = random.Random(seed)
    parent = [-1]*n
    state = [False]*n
    zeros = [0]*n
    total = maximum = 0
    histogram = Counter()
    for x in range(1,n):
        v = parent[x] = rng.randrange(x)
        zeros[v] += 1
        u = k = 0
        while True:
            u += 1
            new = zeros[v] > 0
            if new == state[v]:
                break
            state[v] = new
            k += 1
            if v == 0:
                break
            p = parent[v]
            zeros[p] += -1 if new else 1
            v = p
        total += u
        maximum = max(maximum,u)
        if x >= n//2:
            histogram[k] += 1
    return {'size':n,'seed':seed,'S_over_N':total/n,'max_U':maximum,
            'late_K_histogram':dict(histogram)}

simulations = [simulate(100000,seed) for seed in range(12)]
mean_ratio = sum(x['S_over_N'] for x in simulations)/len(simulations)
print('simulation mean S/N', mean_ratio, flush=True)
report={'exhaustive_trees':checked,'records':records,
        'exact_checks':['history identities','root history mean',
                        'weighted fixed-size covariance','fringe path mean and variance',
                        'pathwise longest-update lower bound','finite PGFs'],
        'simulations':simulations}
Path(__file__).with_name('verification_results.json').write_text(
    json.dumps(report,indent=2),encoding='utf-8')
print('ALL CHECKS PASSED; simulated data are supporting checks, not proofs.',flush=True)
