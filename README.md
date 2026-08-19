# LitRaPID

LitRaPID is a literature- and patent-driven virtual macrocyclic peptide discovery platform. It structures public evidence, represents non-natural macrocycles, builds target evidence packs, generates constrained virtual libraries, and ranks candidates with explicit uncertainty and patent-similarity alerts.

## Repository contents

- `docs/平台设计说明书.md` — complete Chinese platform specification.
- `data/core_data_model.csv` — proposed knowledge-base schema.
- `data/candidate_output_template.csv` — candidate output contract.
- `slides/LitRaPID_Overview_EN.pptx` — English summary presentation.

## Product boundary

LitRaPID prioritizes computational candidates. Without experimental validation, it does not establish binding, selectivity, permeability, safety, or pharmacological activity. Patent outputs are alerts for professional review, not freedom-to-operate opinions.

## Evidence sources

- Public mRNA-display studies and supplementary datasets
- Public patent publications and claim text
- ChEMBL target, assay, molecule, and activity records
- Protein sequence and structural evidence linked to normalized target IDs

## Status

Concept and data-model design. The next implementation milestone is a curated evidence MVP covering 5–10 targets and 500–2,000 reviewed structure–activity records.
