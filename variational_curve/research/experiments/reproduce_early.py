"""Rebuild the early study's numerical tables and publication plots from the preserved generator."""
import argparse
import json
from pathlib import Path
import sys
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code" / "early_studies" / "scripts"))
import finalize_numerical_study as early


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "code" / "results" / "early_surface")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    with threadpool_limits(limits=1):
        result = early.geometry(args.output)
    (args.output / "reproduced.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf8")
    print(args.output)
if __name__=="__main__":main()
