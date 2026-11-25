"""Model families.  Each module defines

    NAME        family name
    LEVELS      difficulty levels it can produce
    FEATURES    {flag: description} -- optional features that ``generate(..., features={flag: bool})``
                can force on or off
    build(draw, level, unit_system, vertical, features) -> Builder

Every family stays inside the envelope that both foreign solvers of ``benchfea.crosscheck`` can
represent (no shear deformation, no notional loads, no stiffness reduction, no rigid offsets in a
model with second-order combinations, diaphragm slaves in the master's plane, no node on the axis
of a member it is not connected to).  Three features leave that envelope and are therefore OPT-IN
flags that are never drawn by default: ``tensionOnly`` bracing, a ``buckling`` analysis in the file
(no foreign solver reproduces either) and ``secondOrder`` in ``diaphragm_frame`` (PyNite's emulation
of a diaphragm does not follow in second order); the metadata says so (``tierC_expected``).
"""
import importlib

FAMILIES = {}
for _name in ("continuous_beam", "portal_frame", "multistorey_2d", "multistorey_3d", "braced_frame", "pipe_rack",
              "stair_frame", "diaphragm_frame"):
    try:
        _m = importlib.import_module("." + _name, __name__)
    except ModuleNotFoundError:
        continue
    FAMILIES[_m.NAME] = _m

__all__ = ["FAMILIES"]
