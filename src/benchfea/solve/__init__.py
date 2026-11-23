"""Solver: ``benchfea.solve.analyze(model) -> Results``.

Modules: ``beam`` (exact integration of the member equations), ``element`` (member stiffness,
releases, offsets, fixed-end forces, recovery), ``loads`` (load sets), ``system`` (assembly,
constraints, supports, factorisation), ``stability`` (mechanism diagnosis), ``invariants``
(self-checks), ``analysis`` (driver and analysis-type dispatch); ``beamcolumn`` (second-order member
piece), ``refine`` (internal subdivision), ``nonlinear`` (second-order iteration, tension/compression-
only members, stiffness reduction, notional loads), ``buckling`` (linear buckling).
Theory: docs/THEORY.md.
"""
from .analysis import ANALYSIS_TYPES, LinearAnalysis, analyze, station_layout, unsupported_features
from .beam import SectionFlex, SpanLoads, integrate
from .element import MemberElement, StiffnessFactors, build_element, cantilever_flexibility
from .errors import MechanismDiagnosis, MechanismMode, ModelRefused, NodeParticipation
from .invariants import model_checks, reciprocity_check, solution_checks, superposition_check
from .loads import LoadSet, applied_resultant, case_loads, combination_loads, combine
from .beamcolumn import BeamColumn, consistent_geometric_stiffness
from .buckling import buckling_analysis
from .nonlinear import CombinationRefused, NonlinearSolver, notional_forces, tau_b
from .refine import Mesh, refine, split_loads
from .system import Factor, Solution, Structure
from ..problems import ModelUnstable

__all__ = ["analyze", "LinearAnalysis", "ANALYSIS_TYPES", "station_layout", "unsupported_features",
           "SectionFlex", "SpanLoads", "integrate", "MemberElement", "StiffnessFactors", "build_element",
           "cantilever_flexibility", "MechanismDiagnosis", "MechanismMode", "ModelRefused", "ModelUnstable",
           "NodeParticipation", "model_checks", "solution_checks", "superposition_check", "reciprocity_check",
           "LoadSet", "applied_resultant", "case_loads", "combination_loads", "combine", "Factor", "Solution",
           "Structure", "BeamColumn", "consistent_geometric_stiffness", "buckling_analysis", "CombinationRefused",
           "NonlinearSolver", "notional_forces", "tau_b", "Mesh", "refine", "split_loads"]
