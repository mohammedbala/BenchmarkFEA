"""Curve ball: second-order effects switched off where they govern."""
from __future__ import annotations

from .registry import Candidate, register

CATEGORY = "second_order"


def _second_order_ids(doc) -> list:
    default = doc.get("analysis", {}).get("type", "firstOrder")
    return [cid for cid, c in doc["combinations"].items() if c.get("analysisType", default) == "secondOrder"]


@register(
    name="second_order_off", category=CATEGORY, title="second-order effects switched off",
    changes="Switches every second-order combination to first order (`/combinations/<id>/analysisType` removed "
            "or set to firstOrder; `/analysis/type` if the file sets it there).",
    why="Under axial compression a frame that sways (P-Delta) and members that bow (P-delta) pick up additional "
        "moments and displacements that a first-order analysis does not see. Where the amplification is material "
        "the first-order results are unconservative: smaller drift, smaller column and beam moments. AISC 360-22 "
        "C1 requires second-order effects to be considered in the required strengths.",
    moves="Drift and moments of the sway-resisting members (down by the amplification), design ratios, the "
          "design status (`conditional`: strength combinations analysed first order).",
    effect_classes=("changes-analysis",),
    precondition="the model has second-order combinations whose results differ materially from first order "
                 "(measured: analysis change of at least 5 %; a change of a design ratio alone does not count)",
    criteria=("analysisChange",))
def second_order_off(ctx):
    from ..difficulty import amplification
    ids = _second_order_ids(ctx.doc)
    if not ids:
        return []
    amps = amplification(ctx.doc, ctx.results)

    def edit(d):
        analysis = d.get("analysis", {})
        inherited = analysis.get("type", "firstOrder") == "secondOrder"
        for cid in ids:
            c = d["combinations"][cid]
            if "analysisType" in c:
                if inherited:
                    c["analysisType"] = "firstOrder"
                else:
                    del c["analysisType"]
        if inherited:
            analysis["type"] = "firstOrder"
    return [Candidate(f"{len(ids)} second-order combination(s) made first order: {', '.join(ids)}", edit,
                      {"combinations": ids},
                      {"baseAmplification": {k: round(float(v), 4) for k, v in amps.items()},
                       "amplificationDefinition": "largest horizontal nodal translation, second order / first order "
                                                  "(benchfea.bench.amplification)"})]
