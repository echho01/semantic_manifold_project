"""Build the self-contained IEEE manuscript and publish only the verified PDF."""
import argparse
from pathlib import Path
import os
import shutil
import subprocess
import tempfile


def run(command, work):
    proc = subprocess.run(command, cwd=work, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if proc.returncode:
        details = proc.stdout.decode(errors="replace").splitlines()
        raise SystemExit("\n".join(details[-45:]))


def main():
    argparse.ArgumentParser(description=__doc__).parse_args()
    project = Path(__file__).resolve().parent
    output = project / "build"
    output.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="semantic-ieee-") as name:
        work = Path(name)
        for filename in ("main.tex", "ref.bib", "IEEEtran.cls"):
            shutil.copyfile(project / filename, work / filename)
        for filename in ("IEEEabrv.bib", "IEEEtran.bst"):
            if (project / filename).is_file():
                shutil.copyfile(project / filename, work / filename)
        if (project / "figures").is_dir():
            shutil.copytree(project / "figures", work / "figures")
        latex = ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", "-file-line-error", "main.tex"]
        run(latex, work)
        run(["bibtex", "main"], work)
        run(latex, work)
        run(latex, work)
        log = (work / "main.log").read_text("utf-8", errors="replace")
        bib = (work / "main.blg").read_text("utf-8", errors="replace")
        if "undefined references" in log or "undefined citations" in log or "Warning--I didn't find" in bib:
            raise SystemExit("Unresolved citation or reference in the manuscript.")
        for line in log.splitlines():
            if "Overfull " in line:
                print(line)
        # Put the temporary file on the destination volume for atomic publication.
        with tempfile.NamedTemporaryFile(dir=output, suffix=".tmp", delete=False) as handle:
            staging = Path(handle.name)
        try:
            shutil.copyfile(work / "main.pdf", staging)
            os.replace(staging, output / "main.pdf")
        finally:
            staging.unlink(missing_ok=True)
    print(output / "main.pdf")


if __name__ == "__main__":
    main()
