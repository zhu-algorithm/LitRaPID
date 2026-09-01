"""Display-ready sequence export and iterative feedback for LitRaPID.

This module prepares computational candidates for Aga2 yeast display,
cysteine-free GPI-anchored yeast macrocycle display, or M13 pIII phage display.
It does not design wet-lab protocols and it does not claim experimental activity.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
from dataclasses import dataclass
from pathlib import Path


AA = set("ACDEFGHIKLMNPQRSTVWY")
HYDROPHOBIC = set("AVILMFWY")
BASIC = set("KRH")
ACIDIC = set("DE")

CODONS = {
    "yeast_aga2": {
        "A":"GCT","C":"TGT","D":"GAT","E":"GAA","F":"TTT","G":"GGT","H":"CAT",
        "I":"ATT","K":"AAA","L":"TTG","M":"ATG","N":"AAT","P":"CCT","Q":"CAA",
        "R":"AGA","S":"TCT","T":"ACT","V":"GTT","W":"TGG","Y":"TAT",
    },
    "m13_p3": {
        "A":"GCG","C":"TGC","D":"GAT","E":"GAA","F":"TTC","G":"GGC","H":"CAC",
        "I":"ATC","K":"AAA","L":"CTG","M":"ATG","N":"AAC","P":"CCG","Q":"CAG",
        "R":"CGT","S":"AGC","T":"ACC","V":"GTG","W":"TGG","Y":"TAC",
    },
}
CODONS["yeast_gpi_macrocycle"] = CODONS["yeast_aga2"]


@dataclass(frozen=True)
class Backend:
    name: str
    min_length: int
    max_length: int
    max_hydrophobic_fraction: float
    max_abs_charge_fraction: float


BACKENDS = {
    "yeast_aga2": Backend("yeast_aga2", 6, 40, 0.50, 0.45),
    # Literature-derived window spanning CX7C/CX9C and CXmCXnC libraries.
    "yeast_gpi_macrocycle": Backend("yeast_gpi_macrocycle", 9, 15, 0.50, 0.45),
    "m13_p3": Backend("m13_p3", 6, 24, 0.45, 0.35),
}

SUBSTITUTIONS = {
    "A": "GSV", "C": "S", "D": "EN", "E": "DQ", "F": "YW", "G": "AS",
    "H": "NQK", "I": "VLM", "K": "RQH", "L": "IVM", "M": "LIV", "N": "QDH",
    "P": "AG", "Q": "NEH", "R": "KQH", "S": "TAG", "T": "SAV", "V": "ILTA",
    "W": "FY", "Y": "FW",
}


def as_float(value: str | None, default: float = 0.0) -> float:
    try:
        return float(value) if value not in (None, "") else default
    except ValueError:
        return default


def logistic(x: float) -> float:
    if x >= 0:
        z = math.exp(-x)
        return 1.0 / (1.0 + z)
    z = math.exp(x)
    return z / (1.0 + z)


def logit(p: float) -> float:
    p = min(1 - 1e-6, max(1e-6, p))
    return math.log(p / (1 - p))


def n_glyco_motifs(sequence: str) -> int:
    return sum(
        sequence[i] == "N" and sequence[i + 1] != "P" and sequence[i + 2] in "ST"
        for i in range(max(0, len(sequence) - 2))
    )


def macrocycle_topology(sequence: str) -> tuple[str, float, list[int]]:
    """Classify a genetically encoded disulfide topology from cysteine positions."""
    positions = [i + 1 for i, aa in enumerate(sequence) if aa == "C"]
    cysteines = len(positions)
    if cysteines == 2:
        return "one_ring", 0.90, positions
    if cysteines == 3:
        return "two_ring_with_random_fourth_cys_possible", 0.72, positions
    if cysteines == 4:
        return "two_disulfide_two_ring", 0.88, positions
    if cysteines > 4 and cysteines % 2 == 0:
        return "multi_disulfide_review", 0.55, positions
    return "topology_review", 0.25, positions


def displayability(sequence: str, backend: Backend) -> tuple[float, list[str], str]:
    seq = sequence.strip().upper()
    reasons: list[str] = []
    if not seq or any(residue not in AA for residue in seq):
        return 0.0, ["contains_noncanonical_or_invalid_residue"], "SYNTHESIS_ONLY"

    score = 1.0
    if len(seq) < backend.min_length or len(seq) > backend.max_length:
        reasons.append("outside_configured_length_window")
        score -= 0.35

    hydrophobic = sum(a in HYDROPHOBIC for a in seq) / len(seq)
    if hydrophobic > backend.max_hydrophobic_fraction:
        reasons.append("high_hydrophobic_fraction")
        score -= min(0.4, hydrophobic - backend.max_hydrophobic_fraction + 0.15)

    charge_fraction = abs(sum(a in BASIC for a in seq) - sum(a in ACIDIC for a in seq)) / len(seq)
    if charge_fraction > backend.max_abs_charge_fraction:
        reasons.append("high_absolute_charge_fraction")
        score -= min(0.35, charge_fraction - backend.max_abs_charge_fraction + 0.15)

    if backend.name == "m13_p3" and seq.count("C") % 2:
        reasons.append("odd_cysteine_count_for_p3")
        score -= 0.35

    if backend.name.startswith("yeast_") and n_glyco_motifs(seq):
        reasons.append("potential_n_glycosylation_motif")
        score -= min(0.3, 0.1 * n_glyco_motifs(seq))

    if backend.name == "yeast_gpi_macrocycle":
        topology, topology_score, _ = macrocycle_topology(seq)
        if seq.count("C") < 2:
            reasons.append("gpi_macrocycle_requires_at_least_two_cysteines")
            score -= 0.55
        if topology == "topology_review":
            reasons.append("disulfide_topology_requires_review")
        score = 0.80 * score + 0.20 * topology_score

    status = "DIRECT_DISPLAY" if score >= 0.55 else "DISPLAY_REVIEW"
    if backend.name == "yeast_gpi_macrocycle" and (
        "outside_configured_length_window" in reasons
        or "gpi_macrocycle_requires_at_least_two_cysteines" in reasons
        or "disulfide_topology_requires_review" in reasons
    ):
        status = "DISPLAY_REVIEW"
    return max(0.0, min(1.0, score)), reasons, status


def backtranslate(sequence: str, backend: str) -> str:
    return "".join(CODONS[backend][aa] for aa in sequence)


def barcode(candidate_id: str, sequence: str, length: int = 12) -> str:
    digest = hashlib.sha256(f"{candidate_id}|{sequence}".encode()).hexdigest()
    bases = "ACGT"
    return "".join(bases[int(ch, 16) % 4] for ch in digest[:length])


def composite_score(row: dict[str, str], display_score: float) -> float:
    affinity = as_float(row.get("predicted_affinity_score"), 0.5)
    selectivity = as_float(row.get("selectivity_score"), 0.5)
    developability = as_float(row.get("developability_score"), 0.5)
    evidence = as_float(row.get("evidence_domain_score"), 0.5)
    uncertainty = as_float(row.get("uncertainty_score"), 0.5)
    return (
        0.30 * affinity + 0.18 * selectivity + 0.15 * developability
        + 0.20 * display_score + 0.10 * evidence + 0.07 * uncertainty
    )


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def export_panel(candidates: Path, backend_name: str, output: Path, limit: int) -> None:
    backend = BACKENDS[backend_name]
    prepared: list[dict[str, object]] = []
    for row in read_csv(candidates):
        seq = (row.get("canonical_sequence") or row.get("peptide_sequence") or "").upper()
        display_score, reasons, status = displayability(seq, backend)
        total = composite_score(row, display_score)
        candidate_id = row.get("candidate_id") or f"CAND-{len(prepared)+1:05d}"
        topology, topology_score, cysteine_positions = macrocycle_topology(seq) if seq and all(a in AA for a in seq) else ("not_applicable", 0.0, [])
        prepared.append({
            "candidate_id": candidate_id,
            "target_id": row.get("target_id", ""),
            "backend": backend_name,
            "peptide_sequence": seq,
            "display_status": status,
            "displayability_score": round(display_score, 6),
            "display_flags": ";".join(reasons),
            "display_anchor": "cysteine_free_GPI" if backend_name == "yeast_gpi_macrocycle" else ("Aga2" if backend_name == "yeast_aga2" else "M13_pIII"),
            "macrocycle_topology": topology,
            "topology_score": round(topology_score, 6),
            "cysteine_positions": ";".join(map(str, cysteine_positions)),
            "selection_plan": "2x_magnetic_bead_enrichment_then_4x_two_color_FACS" if backend_name == "yeast_gpi_macrocycle" else "backend_defined",
            "prior_model_score": round(total, 6),
            "dna_insert": backtranslate(seq, backend_name) if status != "SYNTHESIS_ONLY" else "",
            "barcode": barcode(candidate_id, seq),
            "cluster_id": row.get("cluster_id", ""),
            "control_type": row.get("control_type", "candidate"),
            "round": 0,
        })

    prepared.sort(key=lambda x: (x["display_status"] == "DIRECT_DISPLAY", x["prior_model_score"]), reverse=True)
    if limit > 0:
        prepared = prepared[:limit]
    fields = list(prepared[0]) if prepared else ["candidate_id"]
    write_csv(output, prepared, fields)


def minmax(values: dict[str, float]) -> dict[str, float]:
    if not values:
        return {}
    lo, hi = min(values.values()), max(values.values())
    if math.isclose(lo, hi):
        return {key: 0.5 for key in values}
    return {key: (value - lo) / (hi - lo) for key, value in values.items()}


def hamming_like(a: str, b: str) -> float:
    length = max(len(a), len(b), 1)
    mismatches = sum(x != y for x, y in zip(a, b)) + abs(len(a) - len(b))
    return mismatches / length


def feedback_update(panel_path: Path, results_path: Path, output: Path, limit: int) -> None:
    panel = {row["candidate_id"]: row for row in read_csv(panel_path)}
    results = read_csv(results_path)
    total_input = sum(as_float(r.get("ngs_input_count"), as_float(r.get("input_count"))) for r in results)
    total_selected = sum(as_float(r.get("ngs_selected_count"), as_float(r.get("selected_count"))) for r in results)
    enrich: dict[str, float] = {}
    binding: dict[str, float] = {}
    expression: dict[str, float] = {}
    off_target: dict[str, float] = {}
    stability: dict[str, float] = {}
    function: dict[str, float] = {}
    binding_display_ratio: dict[str, float] = {}
    topology_consistency: dict[str, float] = {}
    for row in results:
        cid = row["candidate_id"]
        input_count = as_float(row.get("ngs_input_count"), as_float(row.get("input_count")))
        selected_count = as_float(row.get("ngs_selected_count"), as_float(row.get("selected_count")))
        inf = (input_count + 0.5) / (total_input + 0.5 * len(results))
        selfreq = (selected_count + 0.5) / (total_selected + 0.5 * len(results))
        enrich[cid] = math.log2(selfreq / inf)
        binding[cid] = as_float(row.get("target_binding_mfi"), as_float(row.get("binding_signal")))
        expression[cid] = as_float(row.get("display_mfi"), as_float(row.get("expression_signal")))
        binding_display_ratio[cid] = binding[cid] / max(expression[cid], 1e-6)
        off_target[cid] = as_float(row.get("off_target_signal"), 0.5)
        stability[cid] = as_float(row.get("stability_signal"), 0.5)
        function[cid] = as_float(row.get("function_signal"), 0.5)
        topology_consistency[cid] = as_float(row.get("topology_consistency"), 0.5)

    e_norm, b_norm, x_norm = minmax(enrich), minmax(binding), minmax(expression)
    o_norm, s_norm, f_norm = minmax(off_target), minmax(stability), minmax(function)
    ratio_norm, t_norm = minmax(binding_display_ratio), minmax(topology_consistency)
    ranked: list[dict[str, object]] = []
    max_round = max((int(as_float(r.get("round"), 0)) for r in results), default=0)
    for cid, old in panel.items():
        prior = as_float(old.get("prior_model_score"), 0.5)
        evidence = (
            0.30 * e_norm.get(cid, 0.5)
            + 0.25 * ratio_norm.get(cid, 0.5)
            + 0.10 * x_norm.get(cid, 0.5)
            + 0.10 * s_norm.get(cid, 0.5)
            + 0.10 * f_norm.get(cid, 0.5)
            + 0.10 * (1.0 - o_norm.get(cid, 0.5))
            + 0.05 * t_norm.get(cid, 0.5)
        )
        posterior = logistic(logit(prior) + 1.2 * (evidence - 0.5))
        updated = dict(old)
        updated.update({
            "round": max_round,
            "log2_enrichment": round(enrich.get(cid, 0.0), 6),
            "binding_signal_norm": round(b_norm.get(cid, 0.5), 6),
            "expression_signal_norm": round(x_norm.get(cid, 0.5), 6),
            "binding_display_ratio": round(binding_display_ratio.get(cid, 0.0), 6),
            "binding_display_ratio_norm": round(ratio_norm.get(cid, 0.5), 6),
            "topology_consistency_norm": round(t_norm.get(cid, 0.5), 6),
            "off_target_signal_norm": round(o_norm.get(cid, 0.5), 6),
            "stability_signal_norm": round(s_norm.get(cid, 0.5), 6),
            "function_signal_norm": round(f_norm.get(cid, 0.5), 6),
            "posterior_score": round(posterior, 6),
        })
        ranked.append(updated)

    ranked.sort(key=lambda x: x["posterior_score"], reverse=True)
    selected: list[dict[str, object]] = []
    for row in ranked:
        seq = str(row["peptide_sequence"])
        if not selected or min(hamming_like(seq, str(chosen["peptide_sequence"])) for chosen in selected) >= 0.20:
            selected.append(row)
        if limit > 0 and len(selected) >= limit:
            break
    if limit <= 0:
        selected = ranked
    fields = list(selected[0]) if selected else ["candidate_id"]
    write_csv(output, selected, fields)


def topology_prior(sequence: str) -> tuple[float, str]:
    """Return a transparent sequence-level topology prior, not a 3D prediction."""
    cysteines = sequence.count("C")
    if cysteines >= 2 and cysteines % 2 == 0:
        return 0.85, "paired-cysteine constrained-loop hypothesis"
    if cysteines == 0:
        return 0.60, "linear or post-display cyclization hypothesis"
    return 0.35, "unpaired-cysteine topology requires review"


def deterministic_variants(sequence: str, parent_id: str, count: int) -> list[tuple[str, str]]:
    """Generate reproducible single substitutions while preserving cysteines."""
    proposals: list[tuple[str, str]] = []
    positions = [i for i, aa in enumerate(sequence) if aa != "C"]
    positions.sort(key=lambda i: hashlib.sha256(f"{parent_id}|{i}".encode()).hexdigest())
    for position in positions:
        source = sequence[position]
        for replacement in SUBSTITUTIONS.get(source, ""):
            variant = sequence[:position] + replacement + sequence[position + 1:]
            proposals.append((variant, f"{source}{position + 1}{replacement}"))
            if len(proposals) >= count:
                return proposals
    return proposals


def design_next_round(
    ranked_path: Path,
    backend_name: str,
    output: Path,
    limit: int,
    parents: int,
    variants_per_parent: int,
) -> None:
    backend = BACKENDS[backend_name]
    ranked = read_csv(ranked_path)
    ranked.sort(key=lambda row: as_float(row.get("posterior_score"), as_float(row.get("prior_model_score"), 0.5)), reverse=True)
    proposals: list[dict[str, object]] = []
    seen: set[str] = set()

    for parent in ranked[:max(1, parents)]:
        parent_id = parent["candidate_id"]
        parent_seq = parent["peptide_sequence"].upper()
        parent_score = as_float(parent.get("posterior_score"), as_float(parent.get("prior_model_score"), 0.5))
        variants = [(parent_seq, "parent-retained")] + deterministic_variants(parent_seq, parent_id, variants_per_parent)
        for index, (sequence, mutation) in enumerate(variants):
            if sequence in seen:
                continue
            seen.add(sequence)
            display_score, flags, status = displayability(sequence, backend)
            topology_score, hypothesis = topology_prior(sequence)
            novelty = hamming_like(sequence, parent_seq)
            acquisition = 0.58 * parent_score + 0.22 * display_score + 0.12 * topology_score + 0.08 * novelty
            next_round = int(as_float(parent.get("round"), 0)) + 1
            cid = f"{parent_id}-R{next_round}-V{index:03d}"
            proposals.append({
                "candidate_id": cid,
                "parent_id": parent_id,
                "target_id": parent.get("target_id", ""),
                "backend": backend_name,
                "peptide_sequence": sequence,
                "mutation_set": mutation,
                "cyclization_topology": parent.get("cyclization_topology", "unspecified"),
                "structure_hypothesis": hypothesis,
                "topology_prior_score": round(topology_score, 6),
                "display_status": status,
                "displayability_score": round(display_score, 6),
                "display_flags": ";".join(flags),
                "parent_posterior_score": round(parent_score, 6),
                "acquisition_score": round(acquisition, 6),
                "dna_insert": backtranslate(sequence, backend_name) if status != "SYNTHESIS_ONLY" else "",
                "barcode": barcode(cid, sequence),
                "round": next_round,
            })

    proposals.sort(key=lambda row: (row["display_status"] == "DIRECT_DISPLAY", row["acquisition_score"]), reverse=True)
    selected: list[dict[str, object]] = []
    for row in proposals:
        seq = str(row["peptide_sequence"])
        if not selected or min(hamming_like(seq, str(chosen["peptide_sequence"])) for chosen in selected) >= 0.10:
            selected.append(row)
        if limit > 0 and len(selected) >= limit:
            break
    fields = list(selected[0]) if selected else ["candidate_id"]
    write_csv(output, selected, fields)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    export = sub.add_parser("export", help="Prepare display-ready sequences")
    export.add_argument("--candidates", required=True, type=Path)
    export.add_argument("--backend", required=True, choices=sorted(BACKENDS))
    export.add_argument("--out", required=True, type=Path)
    export.add_argument("--limit", type=int, default=96)
    feedback = sub.add_parser("feedback", help="Update ranking from display results")
    feedback.add_argument("--panel", required=True, type=Path)
    feedback.add_argument("--results", required=True, type=Path)
    feedback.add_argument("--out", required=True, type=Path)
    feedback.add_argument("--limit", type=int, default=96)
    design = sub.add_parser("design", help="Generate a diverse next-round variant panel")
    design.add_argument("--ranked", required=True, type=Path)
    design.add_argument("--backend", required=True, choices=sorted(BACKENDS))
    design.add_argument("--out", required=True, type=Path)
    design.add_argument("--limit", type=int, default=96)
    design.add_argument("--parents", type=int, default=8)
    design.add_argument("--variants-per-parent", type=int, default=12)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.command == "export":
        export_panel(args.candidates, args.backend, args.out, args.limit)
    elif args.command == "feedback":
        feedback_update(args.panel, args.results, args.out, args.limit)
    else:
        design_next_round(args.ranked, args.backend, args.out, args.limit, args.parents, args.variants_per_parent)


if __name__ == "__main__":
    main()
