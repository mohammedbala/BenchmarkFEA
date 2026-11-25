"""Verification pipeline: decides whether a reference solution may be published as verified.

    from benchfea.verify import verify_model, verify_problem, verify_library
    record = verify_model("frame.json")                # VerificationRecord
    record.status                                       # "verified" | "quarantined" | "refused"
    record.tiers["C"].passed, record.tiers_passed, record.quarantine, record.to_dict()

Tiers (docs/VERIFICATION.md): A closed form / ODE oracle, B published problems, C two independent
solvers, D invariants, E equivalence checks of the modelling features, F coverage of the design checks.
A reference solution is "verified" only if it passed every tier that applies; the record lists them.
"Verified" means exactly that -- agreement of computations at the stated tolerances -- and nothing more.
"""
from .pipeline import software_versions, verify_model
from .policy import POLICY, printed_precision_tolerance
from .problems import (EXTRACTORS, ProblemError, ProblemSpec, find_problems, library_root, load_problem,
                       verify_library, verify_problem)
from .record import (MEANING, QUARANTINED, REFUSED, TIERS, VERIFIED, LibraryReport, ProblemRecord, TargetResult,
                     TierResult, VerificationRecord)
from .tier_c import ForeignCache
from .tier_f import clause_items, coverage_registry, item_coverage

__all__ = ["verify_model", "verify_problem", "verify_library", "VerificationRecord", "TierResult", "ProblemRecord",
           "TargetResult", "LibraryReport", "VERIFIED", "QUARANTINED", "REFUSED", "TIERS", "MEANING", "POLICY",
           "printed_precision_tolerance", "ForeignCache", "load_problem", "find_problems", "library_root",
           "ProblemSpec", "ProblemError", "EXTRACTORS", "software_versions", "coverage_registry", "item_coverage",
           "clause_items"]
