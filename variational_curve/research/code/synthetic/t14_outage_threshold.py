#!/usr/bin/env python3
"""Compatibility entry point for the reproducible Paper 1 experiment suite."""
from pathlib import Path
from threadpoolctl import threadpool_limits
from run_paper1 import run
if __name__ == "__main__":
    with threadpool_limits(limits=1):
        run(Path(__file__).resolve().parents[1] / "results" / "paper1")
