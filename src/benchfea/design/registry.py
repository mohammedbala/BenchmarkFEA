"""Registry of design codes and editions.

An edition is a named tuple of limit-state functions (see
:mod:`benchfea.design.base` for their contract).  Only editions whose text was
supplied and implemented are registered; asking for anything else is refused
(:class:`UnknownDesignCode`) -- there is no fallback to "the nearest" edition.

Adding an edition: create a package next to ``aisc360_22`` with one module per
chapter, each exposing ``LIMIT_STATES``, and call :func:`register` with an
:class:`Edition` from the package's ``__init__``; import the package in
``benchfea/design/__init__.py``.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable


@dataclass(frozen=True)
class Edition:
    code: str                         # the value of the design variable `code`, e.g. "AISC 360-22"
    title: str
    limit_states: tuple               # functions (ctx: DesignContext, demand: Demand) -> list[DesignCheck]
    methods: tuple = ("LRFD", "ASD")
    #: serviceability checks: functions (model, results, design_member_id, ctx) -> list[DesignCheck]
    service: tuple = ()

    def check(self, ctx, demand) -> list:
        """All checks of one design member for one combination's demand."""
        return [c for fn in self.limit_states for c in fn(ctx, demand)]


class UnknownDesignCode(KeyError):
    """The requested design code / edition is not registered."""

    def __init__(self, code):
        self.code = code
        self.registered = registered()
        super().__init__(f"design code {code!r} is not available; registered: {', '.join(self.registered) or 'none'}")

    def __str__(self):
        return self.args[0]


_EDITIONS: dict[str, Edition] = {}


def register(edition: Edition) -> Edition:
    if edition.code in _EDITIONS and _EDITIONS[edition.code] is not edition:
        raise ValueError(f"design code {edition.code!r} is already registered")
    _EDITIONS[edition.code] = edition
    return edition


def get(code: str) -> Edition:
    """The registered edition with exactly this code; raises :class:`UnknownDesignCode`."""
    try:
        return _EDITIONS[code]
    except (KeyError, TypeError):
        raise UnknownDesignCode(code) from None


def registered() -> tuple:
    return tuple(_EDITIONS)


LimitState = Callable[..., list]
