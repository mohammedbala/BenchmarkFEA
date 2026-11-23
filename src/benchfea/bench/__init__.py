"""Benchmark model generators and worked-example support.

    from benchfea.bench import generate, difficulty, families, levels
    g = generate("braced_frame", seed=3, difficulty=3)      # Generated(document, metadata)
    difficulty(g.document)["level"]                         # measured: 3

Families, feature flags and the envelope each family stays in: ``benchfea.bench.families``;
definition of the difficulty levels: docs/DIFFICULTY.md.
"""
from .difficulty import amplification, difficulty, feature_flags, level_of, static_properties
from .families import FAMILIES
from .families._common import CATALOGUE, UNIT_SYSTEMS, Draw, section_definition, stable_seed
from .generate import (BUCKLING_MIN, GENERATOR_VERSION, RATIO_HIGH, RATIO_LOW, RATIO_MAX, Generated,
                       GenerationFailed, families, generate, levels, tier_c_expected)

__all__ = ["generate", "Generated", "GenerationFailed", "families", "levels", "difficulty", "feature_flags",
           "static_properties", "amplification", "level_of", "tier_c_expected", "FAMILIES", "CATALOGUE",
           "UNIT_SYSTEMS", "Draw", "stable_seed", "section_definition", "RATIO_LOW", "RATIO_HIGH", "RATIO_MAX",
           "BUCKLING_MIN", "GENERATOR_VERSION"]
