"""Cross-checks of benchfea against independent solvers (OpenSeesPy, PyNite); docs/CROSSCHECK.md.

    from benchfea.crosscheck import run, compare, cross_check
    r = run("opensees", model, "LC1")              # ForeignResult in benchfea conventions + fidelity record
    report = cross_check(model, "LC1", benchfea_result=results)

Neither foreign package is imported by importing this package; they are optional
(``pip install benchfea[crosscheck]``) and run in a subprocess.
"""
from .capabilities import CAPABILITIES, capability
from .compare import (Comparison, CrossCheckReport, Tolerances, applied_resultant, compare, cross_check,
                      equilibrium, from_benchfea, relative_difference)
from .fidelity import APPROXIMATE, EXACT, UNSUPPORTED, Fidelity, FidelityItem
from .foreign import ForeignResult, Translation, TranslationError
from .runner import SOLVERS, available, convergence_study, run, run_job

__all__ = ["run", "run_job", "available", "convergence_study", "SOLVERS", "compare", "cross_check", "equilibrium", "applied_resultant",
           "from_benchfea", "relative_difference", "Comparison", "CrossCheckReport", "Tolerances", "Fidelity",
           "FidelityItem", "EXACT", "APPROXIMATE", "UNSUPPORTED", "ForeignResult", "Translation",
           "TranslationError", "CAPABILITIES", "capability"]
