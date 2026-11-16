"""One process per declared experiment. plan writes commands; train/validation/test execute that stage."""

import argparse, hashlib, itertools, json, shlex, subprocess, sys
from pathlib import Path


def main():
    here = Path(__file__).resolve().parent
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("stage", choices=("plan", "train", "validation", "test"))
    p.add_argument("--profile", choices=("pilot", "full"), default="pilot")
    p.add_argument("--data", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--device", default="cuda")
    p.add_argument("--download", action="store_true")
    p.add_argument("--spec", type=Path, default=here / "matrix.json")
    a = p.parse_args()
    cfg = json.loads(a.spec.read_text("utf8"))
    spec = cfg["profiles"][a.profile]
    a.out.mkdir(parents=True, exist_ok=True)
    options = list(spec["grid"])
    rows = []
    for values in itertools.product(*(spec["grid"][k] for k in options)):
        axes = dict(zip(options, values))
        for variant in spec.get("variants", [{}]):
            opts = {**spec.get("train", {}), **axes, **variant}
            signature = json.dumps([cfg["entry"], opts], sort_keys=True)
            tag = hashlib.sha256(signature.encode()).hexdigest()[:12]
            out = a.out / a.profile / tag
            ck = out / "codec.pt"
            base = [sys.executable, str(here / cfg["entry"])]

            def argv(obj):
                ret = []
                for k, v in obj.items():
                    ret += ["--" + k, *map(str, v if isinstance(v, list) else [v])]
                return ret

            common = [
                "--device",
                a.device,
                "--data",
                str(a.data.resolve()),
                "--checkpoint",
                str(ck.resolve()),
            ]
            cmds = {"train": base + ["train", *common] + argv(opts)}
            if a.download and cfg["entry"].startswith("semasia"):
                cmds["train"].append("--download")
            for part in ["validation", "test"]:
                sp = (
                    ("val" if part == "validation" else "test")
                    if cfg["entry"] == "topology_jscc.py"
                    else part
                )
                cmds[part] = (
                    base
                    + [
                        "test",
                        *common,
                        "--split",
                        sp,
                        "--output",
                        str((out / (part + ".json")).resolve()),
                    ]
                    + argv(spec.get("eval", {}))
                )
                if a.download and cfg["entry"].startswith("semasia"):
                    cmds[part].append("--download")
            rows.append(dict(id=tag, settings=opts, commands=cmds))
    manifest = a.out / (a.profile + "_matrix.json")
    manifest.write_text(json.dumps(rows, indent=2) + "\n", "utf8")
    print(f"{len(rows)} configurations; manifest: {manifest}", flush=True)
    if a.stage != "plan":
        for row in rows:
            cmd = row["commands"][a.stage]
            print(shlex.join(cmd), flush=True)
            subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
