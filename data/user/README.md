# data/user — files you supply yourself (never committed)

This directory is for data that this project may not distribute. Everything in it except this
README is ignored by Git (`.gitignore`), so a file placed here is not committed by accident.

## The AISC Shapes Database

Needed only if an input file names a section as a standard shape:

```json
"sections": {"col": {"type": "standard", "shape": "W14X90"}}
```

Without the database everything else works: sections given by dimensions (`I`, `box`, `pipe`,
`rectangle`) or by properties need no file. A shape name is then refused with a message that says
what to do; it is never guessed.

**The database is published by AISC. It is not in this repository, and this project gives you no
right to it. Check its licence terms yourself before you copy it or pass it on to anyone.**

### Where to put it

Either of:

* `data/user/aisc-shapes-database-v16.0.xlsx` or `data/user/aisc-shapes-database-v16.0.csv`, in the
  folder you run `benchfea` from or in the repository;
* any path, named by the environment variable `BENCHFEA_SHAPES_DB`.

Reading an `.xlsx` file needs the optional package `openpyxl`: `pip install -e ".[shapes]"`.

### Expected format

The reader was written without access to the database (its download page could not be opened), so
the layout below is an **assumption**. If your file differs, use the two overrides at the end; do
not edit the file.

* The first row of the first worksheet (or of the CSV file) holds the column names.
* Columns read (the first occurrence of a name is used, and is taken to be the US customary
  value):

  | Column | Used as |
  |---|---|
  | `AISC_Manual_Label` or `EDI_Std_Nomenclature` | the shape name, matched exactly (for example `W14X90`, `HSS6X4X1/4`) |
  | `Type` | the kind of shape: `W`, `M`, `S`, `HP` are I-shapes; `HSS`, `PIPE` are hollow sections |
  | `A` | area |
  | `d`, `bf`, `tw`, `tf` | dimensions of an I-shape |
  | `Ht`, `B`, `tdes` | dimensions of a rectangular hollow section |
  | `OD`, `tdes` | dimensions of a round hollow section or pipe |
  | `kdes` | optional; used by the design checks for the web depth of a rolled I-shape |
  | `Ix`, `Zx`, `Sx`, `Iy`, `Zy`, `Sy`, `J`, `Cw` | section properties |

* Lengths are in inches. The values are converted to the units of the input file.

### If your file is laid out differently

* `BENCHFEA_SHAPES_COLUMNS` = path of a JSON file that maps the names above to the column names of
  your file, for example `{"label": ["Shape"], "Ix": "I_major"}`; the logical names are those of
  `COLUMNS` in `src/benchfea/shapes.py` (`label`, `type`, `A`, `d`, `bf`, …).
* `BENCHFEA_SHAPES_LENGTH_UNIT` = `in` or `mm`.

The design checks cover doubly symmetric I-shapes and rectangular and round hollow sections only
([LIMITS.md, section 3.2](../../LIMITS.md#32-code-and-shapes)); other shapes of the database can be
analysed but are not design-checked.
