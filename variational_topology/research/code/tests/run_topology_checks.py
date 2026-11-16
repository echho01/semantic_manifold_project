#!/usr/bin/env python
"""Run independent numerical checks and regression tests; save audit evidence."""
import argparse
import contextlib
import hashlib
import io
import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / 'code/topology'), str(Path(__file__).parent)]
import selftest_tda


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=ROOT / 'code/results/topology/validation')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    log = io.StringIO()
    with contextlib.redirect_stdout(log):
        independent = selftest_tda.run_checks()
    suite = unittest.defaultTestLoader.discover(str(Path(__file__).parent), pattern='test_topology.py')
    result = unittest.TextTestRunner(stream=log, verbosity=2).run(suite)
    sources = list((ROOT / 'code/topology').glob('*.py')) + list((ROOT / 'code/embeddings').glob('*topology*.py')) + [Path(__file__), ROOT / 'code/tests/test_topology.py', ROOT / 'code/embeddings/data_sources.py', ROOT / 'code/embeddings/run_all.py']
    report = {'passed': bool(independent['passed'] and result.wasSuccessful()),
        'independent_libraries': independent, 'regressions': {'run': result.testsRun,
        'failures': len(result.failures), 'errors': len(result.errors)},
        'source_sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}}
    (args.output / 'results.json').write_text(json.dumps(report, indent=2) + '\n')
    (args.output / 'test_report.txt').write_text(log.getvalue())
    print(log.getvalue(), end='')
    return 0 if report['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
