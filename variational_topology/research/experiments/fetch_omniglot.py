"""Fetch official archives by commit and SHA-256; unpack PNGs only."""

import argparse, hashlib, shutil, tempfile, zipfile
from pathlib import Path
from urllib.request import urlopen

REV = "057f034baf2ecb8530bc5710e5a23584d2a519cc"
HASH = {
    "images_background": "ad41ab679c8b5d90b271ef46be6987cc81211774a695c29dcc5367b2b26ee640",
    "images_evaluation": "1f61a8f3366785b057fc117d9228e78a16e3d976c8953b2a10fcc74cf0609cee",
}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data", type=Path, required=True)
    dest = ap.parse_args().data.resolve()
    dest.mkdir(parents=True, exist_ok=True)
    for folder, checksum in HASH.items():
        url = f"https://raw.githubusercontent.com/brendenlake/omniglot/{REV}/python/{folder}.zip"
        with tempfile.TemporaryDirectory() as td:
            archive = Path(td) / "source.zip"
            with urlopen(url, timeout=90) as resp, archive.open("wb") as out:
                shutil.copyfileobj(resp, out)
            if hashlib.sha256(archive.read_bytes()).hexdigest() != checksum:
                raise ValueError("archive hash mismatch")
            count = 0
            with zipfile.ZipFile(archive) as z:
                for member in z.infolist():
                    p = Path(member.filename)
                    if (
                        p.suffix.lower() != ".png"
                        or len(p.parts) != 4
                        or p.parts[0] != folder
                    ):
                        continue
                    out = (dest / p).resolve()
                    if not out.is_relative_to(dest):
                        raise ValueError("invalid archive path")
                    out.parent.mkdir(parents=True, exist_ok=True)
                    with z.open(member) as src, out.open("wb") as dst:
                        shutil.copyfileobj(src, dst)
                    count += 1
            print(folder, count)


if __name__ == "__main__":
    main()
