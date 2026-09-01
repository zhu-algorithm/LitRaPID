# Display-ready sequence and iterative feedback module

## Purpose

This module connects LitRaPID predictions to one explicitly selected display backend:

```text
in silico candidates
  → backend-specific displayability filters
  → peptide + DNA insert + barcode export
  → experimental measurements
  → normalized enrichment / binding / expression / off-target / stability / function evidence
  → posterior ranking
  → constrained mutation design + diversity-aware next round
  → repeat until a predefined stopping rule is met
```

Supported MVP backends:

- `yeast_aga2`: Aga2-based yeast surface display.
- `yeast_gpi_macrocycle`: cysteine-free GPI-anchor yeast display for genetically encoded disulfide macrocycles, based on Linciano et al. (Nature Communications, 2025).
- `m13_p3`: M13 pIII phage display.

They are deliberately separate configurations. A sequence accepted for one backend is not automatically considered suitable for the other.

## Export a display panel

```powershell
$env:PYTHONPATH="src"
python -m litrapid.display_loop export `
  --candidates examples/candidates.csv `
  --backend yeast_aga2 `
  --out work/yeast_panel.csv `
  --limit 96
```

For phage display, use `--backend m13_p3`.

For disulfide macrocycles, use `--backend yeast_gpi_macrocycle`. The exporter records the inferred cysteine topology, cysteine positions, a transparent topology score, the `cysteine_free_GPI` anchor and a literature-derived campaign template of two magnetic-bead enrichment cycles followed by four two-colour FACS cycles. These values are planning metadata, not a wet-lab protocol.

The export contains the peptide sequence, deterministic DNA insert, barcode, displayability score, flags, prior model score and route:

- `DIRECT_DISPLAY`: passes the configured backend rules.
- `DISPLAY_REVIEW`: genetically encodable but contains one or more display risks.
- `SYNTHESIS_ONLY`: contains noncanonical or invalid residue notation and is not directly exportable.

The codon tables are transparent MVP defaults. Production use should replace them with the exact host/vector policy and cloning boundaries used by the laboratory.

## Feed results back

Results use this schema:

```text
candidate_id,round,input_count,selected_count,binding_signal,expression_signal,off_target_signal,stability_signal,function_signal
```

The feedback reader also accepts literature-aligned aliases: `ngs_input_count`, `ngs_selected_count`, `target_binding_mfi`, `display_mfi`, and `topology_consistency`. When present, it calculates a binding/display ratio so affinity evidence is not confounded by surface-expression level.

Run:

```powershell
python -m litrapid.display_loop feedback `
  --panel work/yeast_panel.csv `
  --results examples/display_results.csv `
  --out work/round_2_panel.csv `
  --limit 96
```

The MVP normalizes NGS enrichment, binding/display ratio, display level, off-target, stability, function and topology-consistency signals, updates the prior score, and greedily enforces a minimum sequence distance so the next round is not composed only of near-duplicates.

## Literature-derived yeast macrocycle route

```text
in silico sequence/topology design
  → DNA construct for cysteine-free GPI anchor
  → one-ring / two-ring yeast macrocycle library
  → 2 magnetic-bead enrichment cycles
  → 4 two-colour FACS cycles (target binding + peptide display)
  → Sanger for selected clones + NGS for population deconvolution
  → abundance, sequence-family and binding/display evidence
  → posterior update and next-round design
```

Reference: https://www.nature.com/articles/s41467-025-60907-x

## Design the next round

```powershell
python -m litrapid.display_loop design `
  --ranked work/round_2_panel.csv `
  --backend yeast_aga2 `
  --out work/round_2_variants.csv `
  --parents 8 `
  --variants-per-parent 12 `
  --limit 96
```

The design step retains top parents and generates deterministic conservative substitutions while preserving cysteines. It combines parent posterior, displayability, a transparent topology prior and novelty into an acquisition score. `structure_hypothesis` is explicitly a sequence-level hypothesis, not a predicted 3D structure.

For genuine structure optimization, attach external structural-model outputs—such as structure confidence, interface score, clash/strain score and ensemble stability—to the candidate table, calibrate them against experimental data, and select a Pareto front rather than a single weighted-score winner.

## Recommended loop and stopping rule

```text
Round 0: in silico generation and display export
Round n: display experiment → feedback update → variant design → next panel
Stop: no meaningful improvement across consecutive rounds, uncertainty is calibrated,
      and leading candidates satisfy binding, selectivity, expression, stability,
      function, structural-confidence and IP-review gates.
```

The platform should report a Pareto-optimal candidate set. A single “best structure” is selected only after the project defines which objectives and trade-offs dominate.

## Production extensions

1. Replace heuristic displayability with backend-specific models trained from observed expression and propagation data.
2. Add biological and technical replicate models instead of point estimates.
3. Add negative-selection and homolog-selectivity channels.
4. Separate clone-level and cluster-level uncertainty.
5. Add structure-model and docking adapters around experimentally supported motifs.
6. Version all backend configurations, codon tables and scoring weights.
7. Preserve a holdout set each round to detect model overfitting.

## Boundary

The module prioritizes sequences and prepares computational DNA inserts. It does not establish binding, specify a wet-lab protocol, or provide an FTO opinion. Final constructs must be reviewed against the actual vector, host, fusion orientation, linker, reading frame and laboratory quality system.
