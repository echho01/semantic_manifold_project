#!/usr/bin/env python3
"""Run continuous RD and held-out quantizer regression checks."""
from pathlib import Path
import subprocess
import sys
if __name__ == "__main__":
    raise SystemExit(subprocess.call([sys.executable, str(Path(__file__).resolve().parents[1] / "tests/test_numerics.py"), "RateTests"]))
