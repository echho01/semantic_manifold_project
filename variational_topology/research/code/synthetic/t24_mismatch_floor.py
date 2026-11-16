"""Compatibility entry point for the current Paper 2 mismatch experiment.

Results, per-trial data and provenance are written to code/results/topology.
See run_topology_suite.py for the precise source and receiver model.
"""
import sys
from run_topology_suite import main

if __name__ == "__main__":
    sys.argv[1:1] = ["--stages", "mismatch"]
    raise SystemExit(main())
