"""Linear solver core: ``benchfea.solve.analyze(model) -> Results``.

Modules: ``beam`` (exact integration of the member equations), ``element`` (member stiffness,
releases, offsets, fixed-end forces, recovery), ``loads`` (load sets), ``system`` (assembly,
constraints, supports, factorisation), ``stability`` (mechanism diagnosis), ``invariants``
(self-checks), ``analysis`` (driver and analysis-type dispatch).  Theory: docs/THEORY.md.
"""
from .analysis import ANALYSIS_TYPES, LinearAnalysis, analyze, station_layout, unsupported_features
from .beam import SectionFlex, SpanLoads, integrate
from .element import MemberElement, StiffnessFactors, build_element, cantilever_flexibility
from .errors import MechanismDiagnosis, MechanismMode, ModelRefused, NodeParticipation
from .invariants import model_checks, reciprocity_check, solution_checks, superposition_check
from .loads import LoadSet, applied_resultant, case_loads, combination_loads, combine
from .system import Factor, Solution, Structure
from ..problems import ModelUnstable

__all__ = ["analyze", "LinearAnalysis", "ANALYSIS_TYPES", "station_layout", "unsupported_features",
           "SectionFlex", "SpanLoads", "integrate", "MemberElement", "StiffnessFactors", "build_element",
           "cantilever_flexibility", "MechanismDiagnosis", "MechanismMode", "ModelRefused", "ModelUnstable",
           "NodeParticipation", "model_checks", "solution_checks", "superposition_check", "reciprocity_check",
           "LoadSet", "applied_resultant", "case_loads", "combination_loads", "combine", "Factor", "Solution",
           "Structure"]
