"""Reading and writing input documents (JSON canonical, YAML equivalent).

``load`` returns a :class:`Document`: the plain data (dicts / lists / str /
int / float / bool / None) plus, when the document came from text, the
line/column of every value keyed by its JSON pointer.  Nothing is validated
here except what must be decided while parsing: syntax errors, duplicate
object keys, non-string keys and non-JSON values (NaN, dates, ...).
"""
from __future__ import annotations

import bisect
import json
import math
import re
from dataclasses import dataclass, field
from json.decoder import scanstring
from pathlib import Path
from typing import Any

import yaml

from .problems import ModelInvalid, Problem, pointer


@dataclass
class Document:
    data: Any
    positions: dict[str, tuple[int, int]] = field(default_factory=dict)
    source: str | None = None
    format: str = "dict"  # "json" | "yaml" | "dict"
    parse_problems: list[Problem] = field(default_factory=list)

    def locate(self, path: str) -> tuple[int | None, int | None]:
        """Line/column of ``path``, falling back to the nearest ancestor."""
        p = path
        while True:
            if p in self.positions:
                return self.positions[p]
            if not p:
                return (None, None)
            p = p.rsplit("/", 1)[0]


# --------------------------------------------------------------------------- JSON
_WS = re.compile(r"[ \t\n\r]*")
_NUM = re.compile(r"-?(?:0|[1-9]\d*)(\.\d+)?([eE][+-]?\d+)?")


class _JsonError(Exception):
    def __init__(self, msg, pos):
        self.msg, self.pos = msg, pos


class _JsonParser:
    """Strict RFC 8259 parser that records the position of every value."""

    def __init__(self, text: str):
        self.t = text
        self.starts = [0] + [m.end() for m in re.finditer(r"\n", text)]
        self.positions: dict[str, tuple[int, int]] = {}
        self.problems: list[Problem] = []

    def lc(self, pos: int) -> tuple[int, int]:
        line = bisect.bisect_right(self.starts, pos)
        return line, pos - self.starts[line - 1] + 1

    def ws(self, i: int) -> int:
        return _WS.match(self.t, i).end()

    def parse(self):
        i = self.ws(0)
        val, i = self.value(i, "")
        i = self.ws(i)
        if i != len(self.t):
            raise _JsonError("unexpected text after the end of the document", i)
        return val

    def value(self, i: int, ptr: str):
        t = self.t
        if i >= len(t):
            raise _JsonError("unexpected end of document", i)
        self.positions[ptr] = self.lc(i)
        c = t[i]
        if c == "{":
            obj: dict[str, Any] = {}
            i = self.ws(i + 1)
            if t[i:i + 1] == "}":
                return obj, i + 1
            while True:
                if t[i:i + 1] != '"':
                    raise _JsonError('expected a property name in double quotes', i)
                kpos = i
                try:
                    key, i = scanstring(t, i + 1, True)
                except json.JSONDecodeError as e:
                    raise _JsonError(e.msg, e.pos) from None
                i = self.ws(i)
                if t[i:i + 1] != ":":
                    raise _JsonError("expected ':' after property name", i)
                i = self.ws(i + 1)
                child = ptr + pointer(key)
                if key in obj:
                    line, col = self.lc(kpos)
                    self.problems.append(Problem(
                        child, f"duplicate key {key!r} in the same object", line=line, col=col,
                        hint="Ids and property names must be unique; rename or remove one of them.",
                        code="id.duplicate"))
                    _, i = self.value(i, child + "~duplicate")
                    self.positions.pop(child + "~duplicate", None)
                else:
                    obj[key], i = self.value(i, child)
                i = self.ws(i)
                if t[i:i + 1] == ",":
                    i = self.ws(i + 1)
                    continue
                if t[i:i + 1] == "}":
                    return obj, i + 1
                raise _JsonError("expected ',' or '}'", i)
        if c == "[":
            arr: list[Any] = []
            i = self.ws(i + 1)
            if t[i:i + 1] == "]":
                return arr, i + 1
            while True:
                v, i = self.value(i, f"{ptr}/{len(arr)}")
                arr.append(v)
                i = self.ws(i)
                if t[i:i + 1] == ",":
                    i = self.ws(i + 1)
                    continue
                if t[i:i + 1] == "]":
                    return arr, i + 1
                raise _JsonError("expected ',' or ']'", i)
        if c == '"':
            try:
                return scanstring(t, i + 1, True)
            except json.JSONDecodeError as e:
                raise _JsonError(e.msg, e.pos) from None
        for lit, val in (("true", True), ("false", False), ("null", None)):
            if t.startswith(lit, i):
                return val, i + len(lit)
        m = _NUM.match(t, i)
        if m and m.end() > i:
            s = m.group(0)
            if m.group(1) or m.group(2):
                f = float(s)
                if math.isinf(f):
                    raise _JsonError(f"number {s} is out of range", i)
                return f, m.end()
            return int(s), m.end()
        raise _JsonError("expected a JSON value (NaN/Infinity and comments are not allowed)", i)


def _loads_json(text: str, source: str | None) -> Document:
    if text.startswith("\ufeff"):          # a UTF-8 byte order mark is not part of the document
        text = text[1:]
    p = _JsonParser(text)
    try:
        data = p.parse()
    except RecursionError:
        raise ModelInvalid([Problem("", "JSON document is nested too deeply", code="parse.json",
                                    hint="An input file is at most a few levels deep; check the brackets.")],
                           source) from None
    except _JsonError as e:
        line, col = p.lc(min(e.pos, len(text)))
        raise ModelInvalid([Problem("", f"JSON syntax error: {e.msg}", line=line, col=col, code="parse.json",
                                    hint="Fix the syntax at this position; the file must be strict JSON "
                                         "(double-quoted names, no trailing commas, no comments).")],
                           source) from None
    return Document(data, p.positions, source, "json", p.problems)


# --------------------------------------------------------------------------- YAML
class _YamlLoader(yaml.SafeLoader):
    pass


# PyYAML (YAML 1.1) reads "1e5" as a string; accept the JSON number grammar too
# so that YAML and JSON documents are equivalent.
_YamlLoader.add_implicit_resolver(
    "tag:yaml.org,2002:float",
    re.compile(r"^[-+]?(?:[0-9]+\.?[0-9]*|\.[0-9]+)[eE][-+]?[0-9]+$"),
    list("-+0123456789."),
)


#: numbers a YAML document may contain: the JSON grammar, plus an optional
#: leading "+" and a bare leading/trailing decimal point
_YAML_NUMBER = re.compile(r"^[-+]?(?:(?:0|[1-9][0-9]*)(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][-+]?[0-9]+)?$")


def _loads_yaml(text: str, source: str | None) -> Document:
    if text.startswith("\ufeff"):
        text = text[1:]
    loader = _YamlLoader(text)
    try:
        root = loader.get_single_node()
    except RecursionError:
        loader.dispose()
        raise ModelInvalid([Problem("", "YAML document is nested too deeply", code="parse.yaml",
                                    hint="An input file is at most a few levels deep; check the indentation.")],
                           source) from None
    except yaml.YAMLError as e:
        mark = getattr(e, "problem_mark", None)
        line, col = (mark.line + 1, mark.column + 1) if mark else (None, None)
        raise ModelInvalid([Problem("", f"YAML syntax error: {getattr(e, 'problem', None) or e}", line=line,
                                    col=col, code="parse.yaml", hint="Fix the YAML syntax at this position.")],
                           source) from None
    finally:
        loader.dispose()
    positions: dict[str, tuple[int, int]] = {}
    problems: list[Problem] = []

    def mark(node):
        return node.start_mark.line + 1, node.start_mark.column + 1

    visited: set[int] = set()

    def build(node, ptr, key_mark=None):
        # an alias (*name) refers to the node of its anchor: locate it at its own key
        positions[ptr] = key_mark if id(node) in visited and key_mark else mark(node)
        visited.add(id(node))
        if isinstance(node, yaml.MappingNode):
            out: dict[str, Any] = {}
            for knode, vnode in node.value:
                line, col = mark(knode)
                if knode.tag == "tag:yaml.org,2002:merge":
                    problems.append(Problem(ptr, "YAML merge keys ('<<') are not supported", line=line, col=col,
                                            hint="Write the properties out explicitly.", code="parse.yaml"))
                    continue
                key = loader.construct_object(knode, deep=True) if isinstance(knode, yaml.ScalarNode) else None
                if not isinstance(key, str):
                    problems.append(Problem(
                        ptr, f"key {knode.value!r} is not a string" if isinstance(knode, yaml.ScalarNode)
                        else "mapping key is not a string", line=line, col=col, code="parse.key-type",
                        hint="Ids and property names must be strings; in YAML put quotes around keys "
                             "that look like numbers, booleans or null (e.g. \"1\": ...)."))
                    continue
                child = ptr + pointer(key)
                if key in out:
                    problems.append(Problem(child, f"duplicate key {key!r} in the same object", line=line, col=col,
                                            hint="Ids and property names must be unique; rename or remove one "
                                                 "of them.", code="id.duplicate"))
                    continue
                out[key] = build(vnode, child, (line, col))
            return out
        if isinstance(node, yaml.SequenceNode):
            return [build(v, f"{ptr}/{k}") for k, v in enumerate(node.value)]
        val = loader.construct_object(node, deep=True)
        if isinstance(val, float) and not math.isfinite(val):
            problems.append(Problem(ptr, "NaN/Infinity are not valid numbers", line=mark(node)[0],
                                    col=mark(node)[1], hint="Use a finite number.", code="parse.value"))
        elif isinstance(val, (int, float)) and not isinstance(val, bool) and node.style is None \
                and not _YAML_NUMBER.match(node.value):
            # YAML 1.1 reads 010 as octal (8), 0x10 as 16, 1_000 as 1000 and 1:30 as 5400
            problems.append(Problem(ptr, f"{node.value!r} is not a plain decimal number (YAML would read it "
                                         f"as {val!r})", line=mark(node)[0], col=mark(node)[1], code="parse.value",
                                    hint="Write numbers as in JSON: decimal digits without leading zeros, "
                                         "underscores, colons or 0x/0o prefixes (e.g. 10, 1000, 0.5, 1e3); "
                                         "quote the value if a string is meant."))
        elif val is not None and not isinstance(val, (str, int, float, bool)):
            problems.append(Problem(ptr, f"value of type {type(val).__name__} has no JSON equivalent",
                                    line=mark(node)[0], col=mark(node)[1], code="parse.value",
                                    hint="Use only strings, numbers, booleans, null, lists and mappings "
                                         "(quote date-like strings)."))
            val = str(val)
        return val

    try:
        data = build(root, "") if root is not None else None
    except RecursionError:
        raise ModelInvalid([Problem("", "YAML document is nested too deeply (or an alias refers to itself)",
                                    code="parse.yaml", hint="An input file is at most a few levels deep.")],
                           source) from None
    except yaml.YAMLError as e:
        raise ModelInvalid([Problem("", f"YAML error: {e}", code="parse.yaml",
                                    hint="Use only plain YAML scalars, lists and mappings.")], source) from None
    return Document(data, positions, source, "yaml", problems)


# --------------------------------------------------------------------------- API
def loads(text: str, format: str = "json", source: str | None = None) -> Document:
    """Parse document text. ``format`` is "json" or "yaml"."""
    if format == "json":
        return _loads_json(text, source)
    if format == "yaml":
        return _loads_yaml(text, source)
    raise ValueError(f"unknown format {format!r}; use 'json' or 'yaml'")


def _plain(data: Any, ptr: str, problems: list[Problem]) -> Any:
    """Copy of in-memory data restricted to what a JSON text can hold: the
    same checks the text parsers make (string keys, finite numbers, JSON types)."""
    if isinstance(data, dict):
        out = {}
        for key, val in data.items():
            if not isinstance(key, str):
                problems.append(Problem(ptr, f"key {key!r} is not a string", code="parse.key-type",
                                        hint="Ids and property names must be strings (e.g. \"1\", not 1)."))
                continue
            out[key] = _plain(val, ptr + pointer(key), problems)
        return out
    if isinstance(data, (list, tuple)):
        return [_plain(v, f"{ptr}/{k}", problems) for k, v in enumerate(data)]
    if data is None or isinstance(data, (str, bool, int)):
        return data
    if isinstance(data, float):                 # includes numpy.float64
        if not math.isfinite(data):
            problems.append(Problem(ptr, "NaN/Infinity are not valid numbers", hint="Use a finite number.",
                                    code="parse.value"))
        return data
    for kind in (int, float):                   # numpy and other number types
        try:
            import numbers
            if isinstance(data, numbers.Integral if kind is int else numbers.Real):
                return _plain(kind(data), ptr, problems)
        except (TypeError, ValueError, OverflowError):
            pass
    problems.append(Problem(ptr, f"value of type {type(data).__name__} has no JSON equivalent", code="parse.value",
                            hint="Use only strings, numbers, booleans, None, lists and dicts."))
    return str(data)


def load(source: "str | Path | dict | Document") -> Document:
    """Load a document from a path (.json, .yaml, .yml), a dict, or pass a
    Document through unchanged.  A dict is copied and held to the same rules
    as a text document (string keys, finite numbers, JSON value types)."""
    if isinstance(source, Document):
        return source
    if isinstance(source, dict):
        problems: list[Problem] = []
        return Document(_plain(source, "", problems), parse_problems=problems)
    path = Path(source)
    ext = path.suffix.lower()
    if ext not in (".json", ".yaml", ".yml"):
        raise ModelInvalid([Problem("", f"unsupported file extension {ext!r}", code="parse.extension",
                                    hint="Use .json (canonical) or .yaml/.yml.")], str(path))
    try:
        text = path.read_text(encoding="utf-8")
    except OSError as e:
        raise ModelInvalid([Problem("", f"cannot read file: {e}", code="parse.io",
                                    hint="Check the path.")], str(path)) from None
    except UnicodeDecodeError as e:
        raise ModelInvalid([Problem("", f"the file is not valid UTF-8 text ({e.reason} at byte {e.start})",
                                    code="parse.encoding", hint="Save the file with UTF-8 encoding.")],
                           str(path)) from None
    return loads(text, "json" if ext == ".json" else "yaml", str(path))


def dumps(data: Any, format: str = "json") -> str:
    """Serialise plain data. The canonical form is JSON: 2-space indent, key
    order preserved, UTF-8, no NaN/Infinity, trailing newline."""
    if format == "json":
        return json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    if format == "yaml":
        return yaml.safe_dump(data, sort_keys=False, default_flow_style=None, allow_unicode=True)
    raise ValueError(f"unknown format {format!r}; use 'json' or 'yaml'")


def dump(data: Any, path: "str | Path", format: str | None = None) -> None:
    """Write plain data to ``path`` (format from the extension unless given)."""
    path = Path(path)
    if format is None:
        format = "json" if path.suffix.lower() == ".json" else "yaml"
    path.write_text(dumps(data, format), encoding="utf-8")
