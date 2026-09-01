# PD-L1 epitope design → LitRaPID-DT validation

## Executable contract

The upstream repository `zhu-algorithm/PD-L1-cyclized-peptide` selects one provenance-bearing epitope profile and exports `litrapid.pdl1-report.v1`. LitRaPID consumes it through `POST /api/pdl1-validation`.

1. Generate and rank PD-L1 candidates under `pdl1_antibody_pd1_facing_v_domain` or `pdl1_patent_mapped_contacts`.
2. Preserve the source sequence, conceptual cyclization, epitope profile, binding/selectivity proxies and source rank.
3. Translate the candidate to the RaPID chemistry route. A sequence without a downstream cysteine receives a terminal Cys in `rapid_display_sequence`; this is explicitly recorded as a changed chemical entity.
4. Simulate puromycin fusion, FIT translation, thioether cyclization, PD-L1 occupancy, washing, counterselection, RT-PCR and NGS over multiple rounds.
5. Return raw and bias-corrected rankings for orthogonal PD-L1 binding and PD-1/PD-L1 blockade validation.

## Target and epitope safety rules

- A PD-1 profile is rejected by the PD-L1 validation endpoint; targets cannot be mixed in one ranking.
- Epitope compatibility is a triage prior, not proof of binding-site overlap.
- Terminal-Cys conversion requires re-docking/re-scoring and synthesis of the actual thioether macrocycle.
- Publishable validation requires competition SPR/BLI, PD-1/PD-L1 blockade, round-resolved NGS and appropriate no-target/off-target controls.

## Minimal request

```json
{
  "pdl1_report": {"schema": "litrapid.pdl1-report.v1", "pdl1_report": {}},
  "parameters": {"rounds": 6, "target_concentration_nm": 20, "ngs_reads": 20000, "seed": 2026}
}
```
