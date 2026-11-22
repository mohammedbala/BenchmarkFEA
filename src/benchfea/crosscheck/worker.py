"""Worker process of ``benchfea.crosscheck.runner``:  ``python -m benchfea.crosscheck.worker job.json out.json``.

Reads a job (JSON), executes it with the foreign solver named in it and writes the raw results (JSON).
Kept apart from ``runner`` so that running it as ``__main__`` does not re-import a package module.
"""
from __future__ import annotations

import importlib
import json
import sys
from pathlib import Path

MODULES = {"opensees": "benchfea.crosscheck.opensees", "pynite": "benchfea.crosscheck.pynite"}


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 2:
        print("usage: python -m benchfea.crosscheck.worker job.json out.json", file=sys.stderr)
        return 2
    job = json.loads(Path(argv[0]).read_text())
    try:
        out = importlib.import_module(MODULES[job["solver"]]).execute(job)
    except Exception as e:                                   # noqa: BLE001 - any library failure is a result
        out = {"ok": False, "error": f"{type(e).__name__}: {e}"}
    Path(argv[1]).write_text(json.dumps(out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
