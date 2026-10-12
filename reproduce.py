"""One-command reproduction, using only Python's standard library and C++17."""
from pathlib import Path
import argparse
import csv
import math
import shlex
import shutil
import subprocess
import sys
import os

ROOT = Path(__file__).resolve().parent


def run(arguments):
    print('+ ' + shlex.join(str(x) for x in arguments), flush=True)
    subprocess.run([str(x) for x in arguments], cwd=ROOT, check=True)


def compare_csv(actual, reference):
    with actual.open(newline='') as stream:
        got = list(csv.DictReader(stream))
    with reference.open(newline='') as stream:
        expected = list(csv.DictReader(stream))
    assert len(got) == len(expected), actual.name
    for left, right in zip(got, expected):
        assert left.keys() == right.keys()
        for key in left:
            a, b = left[key], right[key]
            if a == b:
                continue
            # Integer counts must agree exactly. Only decimal output admits the
            # small tolerance needed for long-double differences across platforms.
            assert any(mark in a + b for mark in ('.', 'e', 'E')), (actual.name, key, a, b)
            assert math.isclose(float(a), float(b), rel_tol=1e-12, abs_tol=1e-12), (actual.name, key, a, b)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--checks-only', action='store_true',
                        help='Run API/brute-force and Python finite checks, without million-vertex experiments.')
    parser.add_argument('--output-dir', type=Path, default=ROOT / 'results',
                        help='Directory for regenerated simulation outputs (default: results/).')
    args = parser.parse_args()
    build = ROOT / 'build'
    build.mkdir(exist_ok=True)
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    compiler = shlex.split(os.environ.get('CXX', 'g++'))
    flags = ['-std=c++17', '-O2', '-Wall', '-Wextra', '-Wpedantic', '-I', ROOT / 'include']
    for source, binary in [('tests/test_matching.cpp', 'test_matching'),
                           ('examples/demo.cpp', 'demo')]:
        run(compiler + flags + [ROOT / source, '-o', build / binary])
    run([build / 'test_matching'])
    run([build / 'demo'])
    run([sys.executable, ROOT / 'verification/verify_extensions.py'])
    run([sys.executable, ROOT / 'verification/verify_preferential.py',
         '--output', output / 'preferential_exact_checks.json'])
    import json
    with (output / 'preferential_exact_checks.json').open() as stream:
        preferential_actual = json.load(stream)
    with (ROOT / 'data/preferential_exact_checks.json').open() as stream:
        assert preferential_actual == json.load(stream), 'Preferential exact results differ'
    if args.checks_only:
        print('All finite and API checks passed.', flush=True)
        return
    run(compiler + flags + [ROOT / 'verification/matching_experiments.cpp',
                            '-o', build / 'matching_experiments'])
    run([build / 'matching_experiments', output])
    run([sys.executable, ROOT / 'verification/summarize_experiments.py', output])
    for name in ('matching_exact_checks.csv', 'matching_runs.csv', 'matching_late_histogram.csv'):
        compare_csv(output / name, ROOT / 'data' / name)
    shutil.copyfile(output / 'experiment_table.tex', ROOT / 'paper/verification/experiment_table.tex')
    print('All regenerated counts and numerical results agree with the saved reference data.', flush=True)


if __name__ == '__main__':
    main()
