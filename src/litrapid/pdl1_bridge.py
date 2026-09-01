"""Bridge PD-L1 epitope-guided design into LitRaPID-DT mRNA display."""
from __future__ import annotations
from dataclasses import asdict, dataclass
from typing import Any
from .mrna_display import SimulationParameters, simulate_mrna_display
from .digital_twin import infer_latent_fitness

@dataclass(frozen=True)
class BridgePolicy:
    terminal_cys_policy: str = "append_if_absent"
    max_display_length: int = 20
    require_pdl1_profile: bool = True

def _rapid_sequence(sequence: str, policy: BridgePolicy) -> tuple[str, str]:
    seq=sequence.strip().upper()
    if not seq or any(x not in "ACDEFGHIKLMNPQRSTVWY" for x in seq): raise ValueError("invalid canonical peptide sequence")
    if seq.endswith("C"): return seq,"unchanged-terminal-cys"
    if policy.terminal_cys_policy == "append_if_absent": return seq+"C","terminal-cys-appended-for-thioether-route"
    if "C" in seq: return seq,"existing-downstream-cys"
    raise ValueError("RaPID thioether route requires a downstream cysteine")

def adapt_pdl1_report(report: dict[str, Any], policy: BridgePolicy | None=None) -> dict[str, Any]:
    if report.get("schema") == "litrapid.pdl1-report.v1": report = report.get("pdl1_report") or {}
    policy=policy or BridgePolicy(); epitope=report.get("epitope_design") or {}; target=str(epitope.get("target", ""))
    if policy.require_pdl1_profile and target.upper().replace("-","") != "PDL1":
        raise ValueError("PD-L1 validation requires a PD-L1 epitope profile; PD-1 profiles cannot be mixed")
    adapted=[]
    for index,row in enumerate(report.get("candidates") or [],1):
        source=str(row.get("sequence") or ""); rapid,conversion=_rapid_sequence(source,policy)
        if len(rapid)>policy.max_display_length: continue
        binding=row.get("binding") or {}; epi=row.get("epitope") or {}; admet=row.get("admet") or {}
        adapted.append({
            "candidate_id":str(row.get("id") or f"PDL1-{index:04d}"),
            "canonical_sequence":rapid,
            "source_sequence":source,
            "source_cyclization":row.get("cyclization","head-to-tail (conceptual)"),
            "display_cyclization":"N-chloroacetyl initiator to downstream Cys thioether",
            "chemistry_conversion":conversion,
            "epitope_profile_id":epitope.get("profile_id",""),
            "epitope_compatibility_proxy":epi.get("antibody_epitope_compatibility_proxy",0),
            "predicted_affinity_score":binding.get("pd_l1_binding_proxy",0.5),
            "selectivity_score":binding.get("selectivity_score",0),
            "developability_score":admet.get("synthesizability",0.5),
            "nonspecific_score":max(0.0,min(1.0,1-float(binding.get("selectivity_score",0)))),
            "source_priority_score":row.get("priority_score",0),
        })
    if len(adapted)<2: raise ValueError("at least two display-compatible PD-L1 candidates are required")
    return {"target":"PD-L1","epitope_profile":epitope,"bridge_policy":asdict(policy),"candidates":adapted,
            "warning":"A terminal-Cys conversion changes the chemical entity; display hits must be re-scored and synthesized in the actual thioether topology."}

def validate_pdl1_in_silico_display(report: dict[str,Any], params: SimulationParameters|None=None,
                                    policy: BridgePolicy|None=None) -> dict[str,Any]:
    bridge=adapt_pdl1_report(report,policy); simulation=simulate_mrna_display(bridge["candidates"],params)
    latent=infer_latent_fitness(simulation["trajectories"])
    return {"bridge":bridge,"simulation":simulation,"bias_corrected_ranking":latent,
            "validation_status":"in-silico process validation; requires PD-L1 binding/blockade experiments"}
