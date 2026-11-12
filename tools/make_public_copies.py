"""Replace personal paths in committed files with public-safe versions.

The first time a file is processed, the original is copied to private/<same path>
(gitignored) and kept there untouched. The committed copy then has:
  - local repo paths (/home/<user>/.../false_positives/) made relative to the repo root
  - the cluster user directory (/data/<group>/<user>) replaced by $DATA, the
    variable ARC sets to that directory (<DATA> in YAML, which does not expand it)

Rerun after regenerating any listed file (e.g. new AlphaJudge summaries):
    python3 tools/make_public_copies.py
"""
import re
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PRIVATE = ROOT / "private"

FILES = [
    "CLUSTER_SETUP.md",
    "3_af3/cluster/fulldb/submit_af3_fulldb.sh",
    "3_af3/cluster/workflow/false_positives.yaml",
    *sorted(str(p.relative_to(ROOT)) for p in (ROOT / "5_alpha_judge").glob("*_summary*.csv")),
]

LOCAL_REPO = re.compile(r"/home/[^/\s,\"]+/[^\s,\"]*?false_positives/")
CLUSTER_USER = re.compile(r"/data/[^/\s]+/[^/\s]+(?=/)")


def sanitise(path: Path, text: str) -> str:
    text = LOCAL_REPO.sub("", text)
    data = "<DATA>" if path.suffix in (".yaml", ".yml") else "$DATA"
    return CLUSTER_USER.sub(data, text)


def main():
    for rel in FILES:
        src = ROOT / rel
        backup = PRIVATE / rel
        if not backup.exists():
            backup.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, backup)
        # always sanitise from the current file so regenerated outputs are picked up
        text = src.read_text(newline="")
        public = sanitise(src, text)
        if public != text:
            src.write_text(public, newline="")
            print(f"sanitised {rel}")


if __name__ == "__main__":
    main()
