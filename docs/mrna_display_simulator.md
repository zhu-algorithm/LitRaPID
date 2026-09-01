# Literature-informed mRNA display simulator

LitRaPID 0.4 implements a stochastic population model of a RaPID-like selection rather than a static ranker. Each candidate is propagated through:

1. DNA transcription and mRNA–puromycin linker ligation;
2. cell-free translation and covalent mRNA–peptide fusion formation;
3. initiator/Cys-compatible thioether macrocyclization;
4. target occupancy, wash survival, nonspecific carry-over and counterselection;
5. reverse transcription and persistent sequence-specific PCR bias;
6. finite-depth NGS sampling, pool regeneration and the next selection round.

The virtual library size defaults to `10^12`; candidates are represented as weighted populations, so the program does not allocate one trillion objects. Binding is calculated using `occupancy = [target] / ([target] + Kd)`. If `kd_nm` is absent, the affinity score becomes a log-scale prior from 10 µM (score 0) to 1 nM (score 1). All efficiencies are configurable and all random sampling is reproducible by seed.

## Evidence used to define the model

- Roberts & Szostak established covalent puromycin-linked genotype–phenotype selection: <https://doi.org/10.1073/pnas.94.23.12297>
- The streamlined protocol enumerates ligation, translation/fusion, selection and reverse-transcription operations and validates enrichment from 10:1–10,000:1 mixtures: <https://pmc.ncbi.nlm.nih.gov/articles/PMC3666848/>
- A RaPID selection used semirandom DNA, mRNA–puromycin conjugation, FIT translation, N-chloroacetyl initiators, spontaneous thioether closure, counterselection, RT-PCR and seven rounds before sequencing: <https://pmc.ncbi.nlm.nih.gov/articles/PMC8189037/>
- Round-by-round high-throughput sequencing showed that four-round mRNA-display enrichment can be quantitatively tracked by sequence abundance: <https://pmc.ncbi.nlm.nih.gov/articles/PMC4563141/>
- Kinetic mRNA display demonstrated that enrichment also contains PCR, transcription, ligation, translation and fusion-formation biases, motivating explicit process-bias terms: <https://pmc.ncbi.nlm.nih.gov/articles/PMC4834215/>

## Scope and interpretation

This is a mechanistic **selection-process simulator** suitable for protocol sensitivity analysis, synthetic benchmarking, NGS trajectory prediction and candidate prioritization. It does not infer a physically correct Kd from sequence alone, simulate ribosome chemistry atom by atom, or predict a final 3D structure. Experimental Kd/kinetic values can be supplied directly; otherwise the score-to-Kd mapping is a prior that must be recalibrated with measured data.

## API

POST `/api/simulate-mrna` with `candidates` and optional `parameters`. Returned data include round summaries, per-sequence trajectories, final ranking, NGS counts, enrichment, occupancy, fusion/cyclization probabilities and PCR bias.
