# data/user — user-supplied data (never committed)

Put your own copy of the **AISC Shapes Database v16.0** here if you want to
refer to sections by standard shape name (`"type": "standard"`).

* Accepted file names: `aisc-shapes-database-v16.0.xlsx` or `.csv`
  (or set the environment variable `BENCHFEA_SHAPES_DB` to any path).
* The database is published by AISC and is **not redistributable**; this
  repository contains no AISC data. Everything in this directory except this
  README is git-ignored.
* Without the file, sections given by properties or by parametric dimensions
  work normally; only shape-name lookup is refused (with a hint).

See `src/benchfea/shapes.py` for the expected columns.
