"""Produce the paper's table from the saved raw experiments (stdlib only)."""
from pathlib import Path
import csv
import json
import math
import statistics
import sys

root = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path(__file__).resolve().parents[1] / 'data'
with (root / 'matching_runs.csv').open(newline='') as stream:
    rows = list(csv.DictReader(stream))
assert len(rows) == 128
summary = []
for n in (1000, 10000, 100000, 1000000):
    group = [r for r in rows if int(r['N']) == n]
    assert sorted(int(r['seed']) for r in group) == list(range(32))
    s = [float(r['S_over_N']) for r in group]
    w = [int(r['W']) for r in group]
    recourse = [float(r['recourse_per_update']) for r in group]
    ell = math.log(n) / math.log(math.log(n))
    assert all(abs(float(r['W_over_ell']) - int(r['W']) / ell) < 1e-12 for r in group)
    summary.append({'N': n, 'runs': 32,
                    'mean_S_over_N': statistics.mean(s),
                    'sd_S_over_N': statistics.stdev(s),
                    'mean_W': statistics.mean(w),
                    'min_W': min(w), 'max_W': max(w),
                    'mean_W_over_ell': statistics.mean(w) / ell,
                    'mean_recourse': statistics.mean(recourse)})

with (root / 'matching_late_histogram.csv').open(newline='') as stream:
    histogram = {int(r['K']): int(r['count']) for r in csv.DictReader(stream)}
count = sum(histogram.values())
assert count == 32 * 500000
mean = sum(k * value for k, value in histogram.items()) / count
variance = sum(k * k * value for k, value in histogram.items()) / count - mean**2
late = {'updates': count, 'mean_K': mean, 'variance_K': variance,
        'P_K_zero': histogram.get(0, 0) / count,
        'limiting_variance': 1.1926947246463882,
        'limiting_P_K_zero': 0.4036526376768059,
        'histogram': histogram}
with (root / 'matching_exact_checks.csv').open(newline='') as stream:
    exact = [{k: int(v) for k, v in r.items()} for r in csv.DictReader(stream)]
assert sum(r['trees'] for r in exact) + 1 == 46234
for r in exact:
    n, trees = r['old_n'], r['trees']
    assert r['sum_K'] == trees
    assert n * r['sum_U'] == (2 * n - 1) * trees
    assert n * r['sum_2_power_K'] == (3 * n - 1) * trees
    assert n * r['root_flips'] == trees
    assert r['sum_recourse'] <= 3 * trees

report = {'model': 'uniform recursive leaf attachment',
          'runs': 32, 'seeds': list(range(32)),
          'rng': 'std::mt19937_64 with explicit rejection sampling',
          'checkpoints': summary, 'late_updates': late,
          'exact_matching_trees_checked': 46234,
          'raw_exact_checks': exact,
          'note': 'Updates within each run are dependent. SD is across runs, not updates.'}
(root / 'experiment_summary.json').write_text(json.dumps(report, indent=2) + '\n')

table = [r'\begin{table}[htbp]', r'\centering\small',
         r'\caption{32 independent growth runs (seeds 0--31), sampled at shared',
         r'checkpoints. SD is the sample standard deviation across runs;',
         r'ranges give the observed extremes of $W_N$. Logarithms are natural.}',
         r'\label{tab:experiments}', r'\begin{tabular}{rccc}', r'\toprule',
         r'$N$ & $S_N/N$: mean (SD) & $W_N$: mean (range) & Mean $W_N/\ell_N$ \\',
         r'\midrule']
for index, r in enumerate(summary):
    table.append(f"$10^{{{index + 3}}}$ & {r['mean_S_over_N']:.4f} ({r['sd_S_over_N']:.4f}) & "
                 f"{r['mean_W']:.2f} ({r['min_W']}--{r['max_W']}) & {r['mean_W_over_ell']:.2f} " + r'\\')
table.extend([r'\bottomrule', r'\end{tabular}', r'\end{table}', ''])
(root / 'experiment_table.tex').write_text('\n'.join(table))
print(json.dumps({'checkpoints': summary, 'late': late}, indent=2))
