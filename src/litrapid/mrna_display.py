"""Stochastic, literature-informed simulator of iterative mRNA display selection.

The simulator follows molecules through fusion formation, RaPID-like thioether
cyclisation, affinity selection, reverse transcription, PCR and NGS sampling.
It is a population-selection model, not an atomistic binding simulator.
"""

from __future__ import annotations

import hashlib
import math
import random
from dataclasses import asdict, dataclass
from typing import Any


@dataclass(frozen=True)
class SimulationParameters:
    rounds: int = 6
    target_concentration_nm: float = 20.0
    virtual_library_size: int = 10**12
    ngs_reads: int = 20_000
    transcription_efficiency: float = 0.92
    puromycin_ligation_efficiency: float = 0.65
    translation_efficiency: float = 0.72
    fusion_efficiency: float = 0.55
    reverse_transcription_efficiency: float = 0.80
    cyclization_efficiency: float = 0.85
    wash_survival: float = 0.65
    nonspecific_retention: float = 2e-4
    counterselection_strength: float = 0.35
    pcr_bias_sigma: float = 0.12
    mutation_rate: float = 0.0
    seed: int = 2025

    def validate(self) -> None:
        if not 1 <= self.rounds <= 20:
            raise ValueError("rounds must be between 1 and 20")
        if self.target_concentration_nm <= 0 or self.virtual_library_size <= 0:
            raise ValueError("target concentration and library size must be positive")
        if not 100 <= self.ngs_reads <= 2_000_000:
            raise ValueError("ngs_reads must be between 100 and 2,000,000")
        for key, value in asdict(self).items():
            if key.endswith("efficiency") or key in {"wash_survival", "counterselection_strength", "mutation_rate"}:
                if not 0 <= value <= 1:
                    raise ValueError(f"{key} must be between 0 and 1")
        if not 0 <= self.nonspecific_retention <= 1 or self.pcr_bias_sigma < 0:
            raise ValueError("invalid retention or PCR-bias parameter")


def _unit_interval(value: Any, default: float) -> float:
    try:
        return min(1.0, max(0.0, float(value)))
    except (TypeError, ValueError):
        return default


def _stable_uniform(label: str, seed: int) -> float:
    digest = hashlib.sha256(f"{seed}:{label}".encode()).digest()
    return (int.from_bytes(digest[:8], "big") + 0.5) / (2**64)


def _candidate_kd_nm(row: dict[str, Any]) -> float:
    raw = row.get("kd_nm") or row.get("estimated_kd_nm")
    if raw not in (None, ""):
        kd = float(raw)
        if kd <= 0:
            raise ValueError("kd_nm must be positive")
        return kd
    score = _unit_interval(row.get("predicted_affinity_score"), 0.5)
    # Log-scale prior: score 0 -> 10 uM, 0.5 -> 100 nM, 1 -> 1 nM.
    return 10 ** (4.0 - 4.0 * score)


def _systematic_multinomial(weights: list[float], reads: int, rng: random.Random) -> list[int]:
    total = sum(weights)
    if total <= 0:
        raise ValueError("all simulated selection weights are zero")
    cumulative, running = [], 0.0
    for weight in weights:
        running += weight / total
        cumulative.append(running)
    counts = [0] * len(weights)
    offset = rng.random() / reads
    index = 0
    for draw in range(reads):
        point = offset + draw / reads
        while index < len(cumulative) - 1 and point > cumulative[index]:
            index += 1
        counts[index] += 1
    return counts


def simulate_mrna_display(candidates: list[dict[str, Any]], params: SimulationParameters | None = None) -> dict[str, Any]:
    """Run a reproducible multi-round mRNA-display population simulation."""
    params = params or SimulationParameters()
    params.validate()
    if len(candidates) < 2:
        raise ValueError("at least two candidates are required")

    rng = random.Random(params.seed)
    prepared: list[dict[str, Any]] = []
    seen: set[str] = set()
    for index, source in enumerate(candidates, 1):
        sequence = str(source.get("canonical_sequence") or source.get("peptide_sequence") or "").strip().upper()
        candidate_id = str(source.get("candidate_id") or f"CAND-{index:04d}")
        if not sequence or any(aa not in "ACDEFGHIKLMNPQRSTVWY" for aa in sequence):
            raise ValueError(f"{candidate_id}: canonical amino-acid sequence required")
        if candidate_id in seen:
            raise ValueError(f"duplicate candidate_id: {candidate_id}")
        seen.add(candidate_id)
        kd = _candidate_kd_nm(source)
        occupancy = params.target_concentration_nm / (params.target_concentration_nm + kd)
        has_cys = "C" in sequence
        cyclization = params.cyclization_efficiency if has_cys else params.cyclization_efficiency * 0.05
        pcr_u = max(1e-12, _stable_uniform(candidate_id + sequence, params.seed))
        pcr_bias = math.exp(params.pcr_bias_sigma * rng.normalvariate(0, 1) + 0.02 * math.log(pcr_u))
        nonspecific = _unit_interval(source.get("nonspecific_score"), _stable_uniform(sequence, params.seed + 9) * 0.25)
        prepared.append({
            "candidate_id": candidate_id,
            "peptide_sequence": sequence,
            "estimated_kd_nm": kd,
            "binding_occupancy": occupancy,
            "cyclization_probability": cyclization,
            "pcr_bias": pcr_bias,
            "nonspecific_score": nonspecific,
        })

    frequencies = [1.0 / len(prepared)] * len(prepared)
    trajectories: list[dict[str, Any]] = []
    summaries: list[dict[str, Any]] = []
    fusion_probability = (
        params.transcription_efficiency
        * params.puromycin_ligation_efficiency
        * params.translation_efficiency
        * params.fusion_efficiency
    )

    for round_number in range(1, params.rounds + 1):
        selection_weights = []
        raw_probabilities = []
        for freq, item in zip(frequencies, prepared):
            specific = item["binding_occupancy"] * params.wash_survival
            background = params.nonspecific_retention * (1 + item["nonspecific_score"])
            counter_factor = max(0.0, 1 - params.counterselection_strength * item["nonspecific_score"])
            selection_probability = fusion_probability * item["cyclization_probability"] * (specific + background) * counter_factor
            raw_probabilities.append(selection_probability)
            selection_weights.append(freq * selection_probability * params.reverse_transcription_efficiency * item["pcr_bias"])

        total = sum(selection_weights)
        selected = [weight / total for weight in selection_weights]
        read_counts = _systematic_multinomial(selected, params.ngs_reads, rng)
        observed = [(count + 0.5) / (params.ngs_reads + 0.5 * len(prepared)) for count in read_counts]
        # NGS-observed composition seeds the next regenerated pool.
        frequencies = observed
        top_index = max(range(len(prepared)), key=lambda i: frequencies[i])
        entropy = -sum(p * math.log(p, 2) for p in frequencies if p > 0)
        summaries.append({
            "round": round_number,
            "top_candidate": prepared[top_index]["candidate_id"],
            "top_frequency": round(frequencies[top_index], 8),
            "observed_candidates": sum(count > 0 for count in read_counts),
            "shannon_entropy_bits": round(entropy, 6),
            "effective_diversity": round(2**entropy, 3),
        })
        for item, selected_freq, count, sel_prob in zip(prepared, selected, read_counts, raw_probabilities):
            # input frequency is reconstructed below from the previous trajectory/initial pool.
            previous = 1 / len(prepared) if round_number == 1 else next(
                row["ngs_frequency"] for row in reversed(trajectories) if row["candidate_id"] == item["candidate_id"]
            )
            enrichment = math.log2((selected_freq + 1e-15) / (previous + 1e-15))
            trajectories.append({
                **item,
                "round": round_number,
                "input_frequency": round(previous, 10),
                "selected_frequency": round(selected_freq, 10),
                "ngs_reads": count,
                "ngs_frequency": round((count + 0.5) / (params.ngs_reads + 0.5 * len(prepared)), 10),
                "log2_enrichment": round(enrichment, 6),
                "fusion_probability": round(fusion_probability, 6),
                "selection_probability": round(sel_prob, 10),
            })

    final_rows = sorted(
        (row for row in trajectories if row["round"] == params.rounds),
        key=lambda row: row["ngs_reads"], reverse=True,
    )
    return {
        "model": "population-level RaPID-like mRNA display",
        "parameters": asdict(params),
        "round_summaries": summaries,
        "trajectories": trajectories,
        "final_ranking": final_rows,
        "limitations": "Population-selection simulator; it does not calculate atomistic structures or replace biochemical validation.",
    }
