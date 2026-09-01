"""Bias-aware inference and counterfactual protocol search for LitRaPID-DT."""
from __future__ import annotations
import math
from dataclasses import replace
from itertools import product
from typing import Any, Iterable
from .mrna_display import SimulationParameters, simulate_mrna_display

def infer_latent_fitness(trajectories: list[dict[str, Any]], ridge: float = 0.25) -> list[dict[str, Any]]:
    """Estimate bias-corrected selection fitness from multiround trajectories.

    PCR bias is treated as a measured control offset; inverse-count shrinkage
    prevents low-read clones from receiving extreme corrections.
    """
    grouped: dict[str, list[dict[str, Any]]] = {}
    for row in trajectories: grouped.setdefault(str(row["candidate_id"]), []).append(row)
    output=[]
    for cid, rows in grouped.items():
        weighted=weight_sum=0.0
        for row in rows:
            reads=max(0,int(row.get("ngs_reads",0))); weight=reads/(reads+25.0)
            corrected=float(row["log2_enrichment"])-math.log2(max(float(row.get("pcr_bias",1)),1e-12))
            weighted += weight*corrected; weight_sum += weight
        estimate=weighted/(weight_sum+ridge)
        output.append({"candidate_id":cid,"latent_fitness":round(estimate,6),"information_weight":round(weight_sum,4),"rounds":len(rows)})
    return sorted(output,key=lambda x:x["latent_fitness"],reverse=True)

def optimize_protocol(candidates: list[dict[str, Any]], base: SimulationParameters,
                      target_grid: Iterable[float], wash_grid: Iterable[float], counter_grid: Iterable[float]) -> list[dict[str, Any]]:
    """Return Pareto-efficient protocols for affinity recovery and diversity."""
    scored=[]
    for target,wash,counter in product(target_grid,wash_grid,counter_grid):
        params=replace(base,target_concentration_nm=float(target),wash_survival=float(wash),counterselection_strength=float(counter))
        result=simulate_mrna_display(candidates,params); final=result["final_ranking"]; summary=result["round_summaries"][-1]
        best_true=min(candidates,key=lambda r:float(r.get("kd_nm",10**9)))["candidate_id"]
        rank=next(i for i,r in enumerate(final,1) if r["candidate_id"]==best_true)
        scored.append({"target_concentration_nm":target,"wash_survival":wash,"counterselection_strength":counter,
                       "best_binder_reciprocal_rank":round(1/rank,4),"effective_diversity":summary["effective_diversity"]})
    pareto=[]
    for item in scored:
        if not any(other["best_binder_reciprocal_rank"]>=item["best_binder_reciprocal_rank"] and
                   other["effective_diversity"]>=item["effective_diversity"] and
                   (other["best_binder_reciprocal_rank"]>item["best_binder_reciprocal_rank"] or
                    other["effective_diversity"]>item["effective_diversity"]) for other in scored): pareto.append(item)
    return sorted(pareto,key=lambda x:(-x["best_binder_reciprocal_rank"],-x["effective_diversity"]))
