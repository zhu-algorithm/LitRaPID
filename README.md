# LitRaPID

## RaPID-like mRNA display simulation (v0.4)

The local platform now runs a reproducible, multi-round population simulation of mRNA–puromycin fusion formation, FIT/RaPID-like thioether cyclization, target binding and washing, nonspecific retention/counterselection, RT-PCR bias and finite-depth NGS sampling. It reports every sequence trajectory and the final enriched pool; see [the model documentation](docs/mrna_display_simulator.md).

This is a literature-informed selection-process simulator, not an atomistic molecular-dynamics engine. Supply measured `kd_nm` values when available; predicted affinity scores are treated as recalibratable priors.

LitRaPID is a literature- and patent-driven virtual macrocyclic peptide discovery platform. It structures public evidence, represents non-natural macrocycles, builds target evidence packs, generates constrained virtual libraries, and ranks candidates with explicit uncertainty and patent-similarity alerts.

## Repository contents

- `docs/平台设计说明书.md` — complete Chinese platform specification.
- `docs/display_feedback_loop.md` — display-ready sequence export and iterative feedback design.
- `data/core_data_model.csv` — proposed knowledge-base schema.
- `data/candidate_output_template.csv` — candidate output contract.
- `data/display_panel_template.csv` — display-ready panel output contract.
- `data/display_feedback_template.csv` — experimental feedback input contract.
- `src/litrapid/display_loop.py` — runnable Aga2, GPI-macrocycle yeast and M13 pIII display-feedback engine.
- `src/litrapid/platform_server.py` + `web/` — zero-dependency local web platform and JSON API.
- `examples/` — candidate and display-result examples.
- `tests/` — unit tests for filtering, export and feedback updates.
- `slides/LitRaPID_Overview_EN.pptx` — English summary presentation.

## Product boundary

LitRaPID prioritizes computational candidates. Without experimental validation, it does not establish binding, selectivity, permeability, safety, or pharmacological activity. Patent outputs are alerts for professional review, not freedom-to-operate opinions.

## Evidence sources

- Public mRNA-display studies and supplementary datasets
- Public patent publications and claim text
- ChEMBL target, assay, molecule, and activity records
- Protein sequence and structural evidence linked to normalized target IDs

## Run the in silico display platform

```powershell
$env:PYTHONPATH="src"
python -m litrapid.platform_server --port 8765
```

Open `http://127.0.0.1:8765`. The UI supports candidate CSV upload, backend-specific display screening, peptide/DNA/barcode export, NGS/FACS posterior updating, next-round variant design and CSV download. All processing is local.

API endpoints:

- `GET /api/meta`
- `POST /api/export`
- `POST /api/feedback`
- `POST /api/design`

## Status

The repository now includes a runnable display-feedback MVP. The next implementation milestone is a curated evidence base covering 5–10 targets and 500–2,000 reviewed structure–activity records, followed by backend-specific model calibration on real display data.
