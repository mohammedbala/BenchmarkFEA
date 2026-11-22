"""Independent reference for single members and continuous beams.

The governing ordinary differential equations of a straight member are integrated directly as a
multi-point boundary value problem with ``scipy.integrate.solve_bvp`` (collocation).  No stiffness
matrix, shape function or fixed-end force appears anywhere in this package, so its numerical path
shares nothing with ``benchfea.solve``.  See :mod:`benchfea.oracle.beam` for the equations and their
derivation.

    from benchfea.oracle import Beam, Bar
    b = Beam(6000.0, EI=8e12, GAs=None)
    b.support(0.0, v="fixed", r="fixed")
    b.dist(0.0, 6000.0, -10.0)
    s = b.solve()
    s.v(6000.0), s.M(0.0), s.reactions[0.0]
"""
from .beam import Bar, BarSolution, Beam, BeamSolution, OracleError

__all__ = ["Beam", "Bar", "BeamSolution", "BarSolution", "OracleError"]
