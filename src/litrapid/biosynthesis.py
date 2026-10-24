"""Mechanistic biosynthesis of displayed macrocycles: mRNA, phage and bacterial routes.

What this adds to LitRaPID
--------------------------
`mrna_display.py` treats biosynthesis as one scalar: `fusion_probability` is the
product of four constants and is identical for every sequence, and cyclisation
is binary on the presence of a cysteine. `display_loop.py` scores phage and
yeast backends with additive heuristic penalties and has no bacterial route.
Neither module asks the question a wet-lab campaign actually answers first:

    *can this platform physically make this molecule at all?*

This module walks the central dogma one step at a time for each platform —
transcription, linker ligation, tRNA acylation, elongation, fusion, ring
closure for the cell-free route; ligation, transformation, amber readthrough,
Sec export, oxidative folding, proteolysis, virion assembly and propagation for
the phage route; transformation, induction, intein splicing, host toxicity and
sorting throughput for the bacterial routes — and makes each step's survival
probability a function of the sequence where the mechanism says it should be.

Three loss categories are tracked separately because they have different
consequences:

    bulk        multiplies yield; costs molecules, not sequence space
    ceiling     hard cap on distinct variants (transformation, sorting)
    sequence    changes *which* variants remain -- the only category that
                silently becomes a prior in any model trained on display output

Diversity is propagated with the Poisson survival of each variant's copies,
P(variant survives) = 1 - exp(-copies * p), so the simulator distinguishes
"we lost 90% of the material" from "we lost 90% of the library".

Status
------
Mechanistically-motivated parameterisation, not a fitted model. Each functional
form encodes a documented mechanism (SN2 macrocyclisation kinetics, the
positive-inside rule, the intein +1 nucleophile requirement, OmpT dibasic
specificity); constants are order-of-magnitude values chosen to reproduce
reported platform-level library sizes and yields. Ranks across sequences and
across platforms are the intended output. Absolute efficiencies are not
calibrated against any single experiment and must be recalibrated on measured
per-step data before they are quoted as yields.
"""

from __future__ import annotations

import argparse
import json
import math
import random
from dataclasses import asdict, dataclass, field
from typing import Any, Callable, Optional

AA = "ACDEFGHIKLMNPQRSTVWY"

# Kyte-Doolittle hydropathy, J Mol Biol 157:105 (1982)
HYDROPATHY = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}
CHARGE = {"K": 1.0, "R": 1.0, "H": 0.1, "D": -1.0, "E": -1.0}
BASIC = frozenset("KR")

#: Non-canonical residues reachable only through flexizyme-charged tRNAs.
#: code -> (property proxy, N-methylated, D-configured)
NONCANONICAL = {
    "y": ("Y", False, True),    # ClAc-D-Tyr, the RaPID initiator
    "a": ("A", False, True),    # D-Ala
    "B": ("A", False, False),   # Aib
    "m": ("L", True, False),    # N-Me-Leu
    "n": ("A", True, False),    # N-Me-Ala
    "f": ("F", True, False),    # N-Me-Phe
}

#: IUPAC degenerate bases used by the supported designs.
IUPAC = {"N": "ACGT", "K": "GT", "S": "CG", "B": "CGT",
         **{b: b for b in "ACGT"}}

#: The standard genetic code, used to derive both the residue distribution and
#: the stop burden of a degenerate design. Deriving these rather than tabulating
#: them avoids the class of error where, say, NNB is credited with the three
#: stops of NNN when B excludes adenine from the third position and leaves only
#: the amber codon.
CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L",
    "CTA": "L", "CTG": "L", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "TCT": "S", "TCC": "S",
    "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S", "CCT": "P", "CCC": "P",
    "CCA": "P", "CCG": "P", "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGG": "W", "CGT": "R", "CGC": "R", "CGA": "R",
    "CGG": "R", "AGA": "R", "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G",
    "GGG": "G", "TAA": "*", "TAG": "*", "TGA": "*",
}

SUPPORTED_DESIGNS = ("NNK", "NNS", "NNN", "NNB", "trimer")


def _expand_codons(codon_design: str) -> list[str]:
    pools = [IUPAC[c.upper()] for c in codon_design]
    return [b1 + b2 + b3 for b1 in pools[0] for b2 in pools[1] for b3 in pools[2]]


def _stop_frequency(codon_design: str) -> float:
    if codon_design == "trimer":
        return 0.0
    codons = _expand_codons(codon_design)
    return sum(1 for c in codons if CODON_TABLE[c] == "*") / len(codons)


#: Stop codons per codon of each degenerate design, enumerated from the code.
#: NNK and NNS carry the amber codon alone (1/32); NNB likewise carries only
#: amber (1/48) because B excludes adenine; NNN carries all three (3/64).
STOP_FREQUENCY = {design: _stop_frequency(design) for design in SUPPORTED_DESIGNS}

#: Heavy atoms each residue contributes to a backbone macrocycle.
ATOMS_PER_RESIDUE = 3


# ---------------------------------------------------------------------------
# Residue-level descriptors
# ---------------------------------------------------------------------------
def residue_properties(residue: str) -> tuple[str, bool, bool]:
    """(property proxy, is N-methylated, is D-configured) for any residue code."""
    if residue in NONCANONICAL:
        return NONCANONICAL[residue]
    return (residue, False, False)


def hydropathy(sequence: str) -> float:
    """Mean Kyte-Doolittle hydropathy; non-canonicals use their property proxy."""
    values = [HYDROPATHY.get(residue_properties(r)[0], 0.0) for r in sequence]
    return sum(values) / len(values) if values else 0.0


def net_charge(sequence: str) -> float:
    return sum(CHARGE.get(residue_properties(r)[0], 0.0) for r in sequence)


def n_methyl_count(sequence: str) -> int:
    return sum(1 for r in sequence if residue_properties(r)[1])


def stop_free_fraction(codon_design: str, n_randomised: int) -> float:
    """Fraction of members with no in-frame stop in the randomised region.

    NNK carries exactly one stop in 32 codons, so a 12-mer randomised region is
    stop-free only (31/32)**12 = 68.3% of the time. Every platform pays this;
    what differs is the amber codon's fate downstream (reassigned in a PURE
    system lacking RF1, read as Gln in a supE host, or terminating).
    """
    per_codon = STOP_FREQUENCY.get(codon_design, 1 / 32)
    return (1.0 - per_codon) ** n_randomised


# ---------------------------------------------------------------------------
# Ring-closure kinetics
# ---------------------------------------------------------------------------
def ring_atom_count(n_residues: int, linker_atoms: int = 0) -> int:
    return n_residues * ATOMS_PER_RESIDUE + linker_atoms


def effective_molarity(n_ring_atoms: int, optimum: float = 27.0,
                       width: float = 0.55) -> float:
    """Relative effective molarity for ring closure, log-normal in ring size.

    Small rings pay enthalpic strain; large rings pay conformational entropy
    (the N**-3/2 Jacobson-Stockmayer limit). A log-normal in ring atoms
    reproduces both lobes with one expression. The default optimum of ~27 heavy
    atoms corresponds to the 8-9 residue thioether rings typical of RaPID.
    """
    if n_ring_atoms <= 6:
        return 1e-4
    return math.exp(-(math.log(n_ring_atoms / optimum) ** 2) / (2 * width ** 2))


def thioether_closure(sequence: str, k_hydrolysis: float = 0.06,
                      k_intermolecular: float = 0.04) -> dict[str, Any]:
    """ClAc-D-Tyr + Cys thiolate SN2 closure, the RaPID chemistry.

    The N-terminal chloroacetamide is attacked by the first downstream cysteine
    thiolate at pH ~8. Competing fates are hydrolysis of the ClAc group and
    intermolecular alkylation of another chain. A second internal cysteine does
    not block closure but splits material between two ring sizes — a real
    liability of cysteine-containing randomised regions that a binary
    "has cysteine" test cannot represent.
    """
    first = sequence.find("C")
    if first < 0:
        return {"efficiency": 0.0, "ring_atoms": 0, "cysteines": 0,
                "isomer_mixture": False, "note": "no cysteine nucleophile"}
    cysteines = sequence.count("C")
    ring = ring_atom_count(first + 1, linker_atoms=2) + 1
    em = effective_molarity(ring)
    efficiency = em / (em + k_hydrolysis + k_intermolecular)
    if cysteines > 1:
        efficiency *= 0.80
    return {"efficiency": efficiency, "ring_atoms": ring,
            "cysteines": cysteines, "isomer_mixture": cysteines > 1, "note": ""}


def disulfide_closure(sequence: str) -> dict[str, Any]:
    """Oxidative closure by DsbA in the periplasm or PDI in the ER.

    Requires a cysteine pair. An odd count leaves a free thiol that drives
    intermolecular scrambling; three or more give a pairing isomer mixture,
    which is why bicyclic phage libraries fix three cysteines and close them
    chemically after the virion is made rather than relying on oxidation.
    """
    cysteines = sequence.count("C")
    if cysteines < 2:
        return {"efficiency": 0.0, "ring_atoms": 0, "cysteines": cysteines,
                "isomer_mixture": False, "note": "fewer than two cysteines"}
    first, last = sequence.find("C"), sequence.rfind("C")
    ring = ring_atom_count(last - first + 1) - 1
    em = effective_molarity(ring, optimum=24.0, width=0.65)
    efficiency = em / (em + 0.10)
    if cysteines == 3:
        efficiency *= 0.45
    elif cysteines > 3:
        efficiency *= 0.25
    return {"efficiency": efficiency, "ring_atoms": ring,
            "cysteines": cysteines, "isomer_mixture": cysteines > 2, "note": ""}


def intein_splicing(extein: str, intein: str = "Npu") -> dict[str, Any]:
    """SICLOPPS head-to-tail closure by a split DnaE intein.

    The extein's first residue must be a nucleophile (Cys >> Ser > Thr) to
    attack the branched intermediate, and the rate is sensitive to the residue
    immediately following it: proline blocks the acyl shift and beta-branched
    residues slow it. Npu DnaE is markedly more tolerant of that position than
    Ssp. Failed splicing routes to N-terminal cleavage, releasing a linear
    by-product rather than a macrocycle.
    """
    if not extein:
        return {"efficiency": 0.0, "ring_atoms": 0, "note": "empty extein"}
    base = {"Npu": 0.72, "Ssp": 0.38}.get(intein, 0.50)
    nucleophile = {"C": 1.0, "S": 0.55, "T": 0.35}.get(extein[0], 0.02)
    tolerance = 1.0 if intein == "Npu" else 0.55
    following = extein[1] if len(extein) > 1 else "G"
    penalty = {"P": 0.15, "I": 0.60, "V": 0.65, "T": 0.70,
               "W": 0.75}.get(following, 1.0)
    penalty = 1.0 - (1.0 - penalty) * tolerance
    ring = ring_atom_count(len(extein))
    em = effective_molarity(ring, optimum=21.0, width=0.50)
    note = "" if nucleophile > 0.1 else "no +1 nucleophile (requires C/S/T)"
    return {"efficiency": base * nucleophile * penalty * (em / (em + 0.08)),
            "ring_atoms": ring, "note": note}


# ---------------------------------------------------------------------------
# Translation and in-cell trafficking
# ---------------------------------------------------------------------------
def elongation_fidelity(sequence: str, canonical: float = 0.997,
                        n_methyl: float = 0.80, d_residue: float = 0.70,
                        consecutive_n_methyl: float = 0.45) -> float:
    """Fraction of ribosomes reaching full length.

    N-methyl amino acids and D-residues are poor A-site substrates, and
    consecutive N-methyls stall the peptidyl transferase centre far worse than
    isolated ones. This is the constraint that stops permeability-driven
    N-methylation patterns from simply being dialled up: the ring that crosses
    the membrane best is the one the ribosome makes worst.
    """
    factor = 1.0
    previous_n_methyl = False
    for residue in sequence:
        _, is_n_methyl, is_d = residue_properties(residue)
        if is_n_methyl:
            factor *= n_methyl * (consecutive_n_methyl if previous_n_methyl else 1.0)
        elif is_d:
            factor *= d_residue
        else:
            factor *= canonical
        previous_n_methyl = is_n_methyl
    return factor


def sec_export(sequence: str) -> float:
    """Sec/SRP translocation competence of a peptide behind a signal sequence.

    Two documented sequence effects dominate. The positive-inside rule: net
    positive charge in the mature N-terminal region opposes the proton-motive
    force and blocks translocation. Strong hydrophobicity diverts the chain
    into the inner membrane or into inclusion bodies. Prolines slow threading
    through SecYEG.
    """
    charge = net_charge(sequence)
    gravy = hydropathy(sequence)
    z = (1.6
         - 0.55 * max(0.0, charge)
         - 0.30 * max(0.0, gravy - 1.2) ** 2
         - 0.25 * sequence.count("P"))
    return 1.0 / (1.0 + math.exp(-z))


def proteolysis_survival(sequence: str) -> float:
    """Survival of OmpT and DegP. OmpT cleaves between consecutive basic residues."""
    dibasic = sum(1 for i in range(len(sequence) - 1)
                  if sequence[i] in BASIC and sequence[i + 1] in BASIC)
    unfolded_load = max(0.0, hydropathy(sequence) - 1.5)
    return 0.92 * (0.35 ** dibasic) * math.exp(-0.25 * unfolded_load)


def host_toxicity(sequence: str) -> float:
    """Probability a clone is lost during propagation.

    Cationic plus hydrophobic peptides are antimicrobial by construction. In a
    library expressed inside the host this removes a biased slice of sequence
    space before any target is introduced, and it removes it from every round.
    """
    charge = max(0.0, net_charge(sequence))
    gravy = max(0.0, hydropathy(sequence))
    return 1.0 / (1.0 + math.exp(-(-2.6 + 0.70 * charge + 0.55 * gravy)))


def n_glycosylation_motifs(sequence: str) -> int:
    """N-X-S/T sequons, which the yeast secretory pathway will glycosylate."""
    return sum(1 for i in range(max(0, len(sequence) - 2))
               if sequence[i] == "N" and sequence[i + 1] != "P"
               and sequence[i + 2] in "ST")


def glycosylation_survival(sequence: str) -> float:
    """Displayed-but-modified fraction; each sequon occludes part of the ring."""
    return 0.75 ** n_glycosylation_motifs(sequence)


# ---------------------------------------------------------------------------
# Platform description
# ---------------------------------------------------------------------------
SequenceTerm = Callable[[str], float]


@dataclass(frozen=True)
class Step:
    """One physical operation in the route from DNA to a displayed macrocycle."""

    token: str
    name: str
    name_zh: str
    layer: str                                  # DNA|RNA|translation|chemistry|linkage|host|readout
    bulk: float = 1.0                           # sequence-independent survival
    ceiling: Optional[float] = None              # hard cap on distinct variants
    sequence_term: Optional[SequenceTerm] = None  # per-sequence survival
    rewrites: Optional[Callable[[str], str]] = None
    mechanism: str = ""
    consequence: str = ""

    @property
    def category(self) -> str:
        if self.ceiling is not None:
            return "ceiling"
        if self.sequence_term is not None:
            return "sequence"
        return "bulk"


@dataclass(frozen=True)
class Platform:
    key: str
    name: str
    name_zh: str
    linkage: str                 # how genotype and phenotype are physically coupled
    compartment: str
    chemistries: tuple[str, ...]
    noncanonical: bool
    readout: str
    steps: tuple[Step, ...]
    #: Molecules of assembled DNA fed into the route. Synthetic DNA is cheap on
    #: every platform, so this is NOT the library size: on cell-based routes
    #: most of this input dies at the transformation step, and the difference
    #: between the two numbers is the whole point of the funnel.
    dna_input: float
    practical_library: float     # reported upper bound on the resulting library
    typical_ring_residues: tuple[int, int]

    def step(self, token: str) -> Step:
        for candidate in self.steps:
            if candidate.token == token:
                return candidate
        raise KeyError(token)


# --- mRNA display (RaPID) --------------------------------------------------
MRNA_DISPLAY = Platform(
    key="mrna_rapid",
    name="mRNA display (RaPID, puromycin-linked)",
    name_zh="mRNA展示（RaPID，嘌呤霉素连接）",
    linkage="covalent peptide-puromycin-mRNA fusion, converted to cDNA",
    compartment="cell-free reconstituted translation (PURE)",
    chemistries=("thioether", "lactam", "head_to_tail", "disulfide"),
    noncanonical=True,
    readout="RT-PCR of the fused cDNA followed by NGS",
    dna_input=1e13,
    practical_library=1e13,
    typical_ring_residues=(6, 15),
    steps=(
        Step("DNA", "DNA library synthesis and PCR", "DNA文库合成与PCR扩增", "DNA",
             bulk=1.0,
             mechanism="degenerate-codon synthesis; no cellular step, so the "
                       "only limit is how many molecules are in the tube",
             consequence="10^13-10^14 molecules is a routine input scale, "
                         "three to four orders above any cell-based route"),
        Step("TXN", "T7 in vitro transcription", "T7体外转录", "RNA",
             bulk=0.85,
             mechanism="abortive initiation and 3' run-off truncation",
             consequence="bulk yield loss only"),
        Step("LIG", "Puromycin-linker ligation", "嘌呤霉素接头连接", "linkage",
             bulk=0.75,
             mechanism="splint-directed T4 RNA ligase joins the puromycin "
                       "linker to the mRNA 3' end",
             consequence="this step CREATES the genotype-phenotype link, so "
                         "its efficiency multiplies everything downstream"),
        Step("ACYL", "Flexizyme charging of the initiator tRNA",
             "柔性酶装载起始tRNA", "translation",
             bulk=0.50,
             mechanism="eFx/dFx ribozyme acylates tRNA-fMet-CAU with "
                       "ClAc-D-Tyr instead of methionine",
             consequence="the step that buys non-canonical chemistry; it has "
                         "no counterpart in any cell-based platform"),
        Step("TRL", "PURE translation without Met and RF1",
             "PURE体系翻译（缺Met/RF1）", "translation",
             bulk=0.60, sequence_term=elongation_fidelity,
             mechanism="reconstituted translation; omitting release factor 1 "
                       "makes the amber codon reassignable rather than "
                       "terminating",
             consequence="N-methyl and D-residues reduce elongation yield, so "
                         "permeability-oriented backbones cost material here"),
        Step("FUS", "Puromycin A-site entry forms the fusion",
             "嘌呤霉素进入A位形成融合体", "linkage",
             bulk=0.25,
             mechanism="the ribosome stalls at the template 3' end and "
                       "puromycin accepts the nascent chain",
             consequence="the dominant bulk loss of the route; 10-40% of "
                         "templates are routinely converted"),
        Step("CYC", "Spontaneous thioether cyclisation", "自发硫醚环化", "chemistry",
             sequence_term=lambda s: thioether_closure(s)["efficiency"],
             mechanism="intramolecular SN2 of the cysteine thiolate on the "
                       "N-terminal chloroacetamide at pH ~8, competing with "
                       "hydrolysis and intermolecular alkylation",
             consequence="ring size and cysteine count decide the yield; "
                         "sequences without a cysteine are never cyclised"),
        Step("RT", "Reverse transcription", "逆转录", "linkage",
             bulk=0.90,
             mechanism="cDNA synthesis removes mRNA secondary structure that "
                       "would otherwise interfere with selection",
             consequence="bulk yield loss only"),
        Step("PUR", "Affinity purification of fusions", "融合体亲和纯化", "readout",
             bulk=0.30,
             mechanism="oligo-dT or FLAG capture of the fusion product",
             consequence="recovery loss; only a fraction of captured fusions "
                         "are conformationally competent"),
    ),
)

# --- phage display (M13 pIII) ----------------------------------------------
PHAGE_DISPLAY = Platform(
    key="phage_m13",
    name="Phage display (M13, pIII fusion)",
    name_zh="噬菌体展示（M13，pIII融合）",
    linkage="peptide on pIII, encoding phagemid packaged in the same virion",
    compartment="E. coli periplasm, then the virion surface",
    chemistries=("disulfide", "chemical_crosslink"),
    noncanonical=False,
    readout="PCR from eluted phage then NGS; amplified by re-infection",
    dna_input=1e13,
    practical_library=1e10,
    typical_ring_residues=(5, 20),
    steps=(
        Step("DNA", "Insert ligation into the phagemid", "插入片段连入噬菌粒", "DNA",
             bulk=0.35,
             mechanism="restriction and ligation of the randomised cassette",
             consequence="bulk loss, and it sets how much DNA reaches "
                         "electroporation"),
        Step("XFORM", "Electroporation into E. coli TG1",
             "电转化进入大肠杆菌TG1", "host",
             bulk=0.05, ceiling=1e10,
             mechanism="each library member must enter a living cell to exist "
                       "at all; 1-10 x 10^9 transformants is the practical "
                       "range and ~10^10 cfu/ug the best case",
             consequence="THE diversity ceiling of every cell-based platform. "
                         "A 12-mer NNK design encodes 4x10^15 sequences, so "
                         "more than 99.999% of the design is never made"),
        Step("SUP", "Amber readthrough in a supE host",
             "supE宿主中的琥珀密码子抑制", "translation",
             bulk=0.60, rewrites=lambda s: s.replace("*", "Q"),
             mechanism="supE tRNA reads TAG as glutamine at roughly 50-70% "
                       "efficiency; the remainder terminates",
             consequence="the amber codon is not lost, it is silently "
                         "rewritten to Gln — a systematic glutamine "
                         "enrichment with no counterpart in mRNA display"),
        Step("EXP", "Sec/SRP export to the periplasm", "Sec/SRP转运至周质", "protein",
             sequence_term=sec_export,
             mechanism="the positive-inside rule and hydrophobicity-driven "
                       "mistargeting govern translocation",
             consequence="cationic and strongly hydrophobic rings are "
                         "depleted before any target is introduced"),
        Step("FOLD", "DsbA oxidative disulfide closure", "DsbA氧化成二硫键", "chemistry",
             sequence_term=lambda s: disulfide_closure(s)["efficiency"],
             mechanism="periplasmic thiol-disulfide oxidoreductase closes a "
                       "cysteine pair supplied by the randomised region",
             consequence="the only genetically encoded cyclisation available "
                         "here; thioether and N-methyl chemistry are not"),
        Step("PROT", "OmpT and DegP proteolysis", "OmpT/DegP蛋白酶降解", "protein",
             sequence_term=proteolysis_survival,
             mechanism="OmpT cleaves between consecutive basic residues; DegP "
                       "clears unfolded passengers",
             consequence="dibasic motifs are selectively destroyed, a "
                         "sequence-space hole specific to this platform"),
        Step("ASSY", "Virion assembly and pIII incorporation",
             "病毒粒子组装与pIII掺入", "linkage",
             bulk=0.12,
             mechanism="a phagemid plus helper phage gives near-monovalent "
                       "display; the Poisson fraction of virions carrying at "
                       "least one fusion copy is around 10%",
             consequence="monovalency is what makes panning affinity-driven "
                         "rather than avidity-driven, and it is paid for in "
                         "yield"),
        Step("PROP", "Propagation and growth-rate bias",
             "增殖过程中的生长速率偏倚", "host",
             sequence_term=lambda s: 1.0 - host_toxicity(s),
             mechanism="clones whose product burdens the host are outgrown "
                       "during every amplification cycle",
             consequence="compounds over rounds: a mildly toxic binder can be "
                         "lost even while it is being enriched"),
    ),
)

# --- bacterial display, intracellular macrocycles (SICLOPPS) ---------------
BACTERIAL_SICLOPPS = Platform(
    key="bacteria_siclopps",
    name="Bacterial SICLOPPS (split-intein head-to-tail cyclisation)",
    name_zh="细菌SICLOPPS（断裂内含肽首尾环化）",
    linkage="macrocycle and its plasmid co-localised inside one cell",
    compartment="E. coli cytoplasm (reducing)",
    chemistries=("head_to_tail",),
    noncanonical=False,
    readout="reverse two-hybrid survival or FACS; plasmid recovered from survivors",
    dna_input=1e12,
    practical_library=1e9,
    typical_ring_residues=(5, 8),
    steps=(
        Step("DNA", "Cassette ligation into the SICLOPPS vector",
             "盒式序列连入SICLOPPS载体", "DNA",
             bulk=0.30,
             mechanism="the I_C-extein-I_N arrangement fixes the extein's "
                       "first residue as cysteine because splicing needs a "
                       "nucleophile there",
             consequence="one of the ring positions is spent on chemistry "
                         "rather than on recognition"),
        Step("XFORM", "Transformation into the reporter strain", "转化入报告菌株",
             "host",
             bulk=0.05, ceiling=1e9,
             mechanism="10^8-10^9 transformants is the reported practical "
                       "SICLOPPS library size",
             consequence="libraries are therefore built at 5-6 randomised "
                         "positions (32^5 = 3.4x10^7), so the ring size that "
                         "fits the platform is set by the ceiling, not by "
                         "chemistry"),
        Step("EXPR", "Arabinose-induced expression", "阿拉伯糖诱导表达", "translation",
             bulk=0.70,
             mechanism="pBAD induction, with cell-to-cell variation in copy "
                       "number and expression level",
             consequence="adds noise to the phenotype of an identical genotype"),
        Step("SPLICE", "Split-intein splicing to the macrolactam",
             "内含肽剪接生成大环内酰胺", "chemistry",
             sequence_term=lambda s: intein_splicing(s, "Npu")["efficiency"],
             mechanism="Npu DnaE splicing; the competing N-terminal cleavage "
                       "pathway releases a linear by-product instead of a ring",
             consequence="the residue following the catalytic cysteine gates "
                         "the yield, with proline nearly abolishing it"),
        Step("TOX", "Host toxicity filter", "宿主毒性筛除", "host",
             sequence_term=lambda s: 1.0 - host_toxicity(s),
             mechanism="an unbiased cyclic-peptide library is substantially "
                       "growth-inhibitory to the host expressing it",
             consequence="a large, charge- and hydrophobicity-biased slice of "
                         "the library is lost before selection begins"),
        Step("SORT", "Reporter readout and sorting throughput",
             "报告系统读出与分选通量", "readout",
             bulk=0.80, ceiling=1e8,
             mechanism="survival selection or FACS caps how many cells can be "
                       "interrogated per round",
             consequence="a second ceiling below the transformation ceiling"),
    ),
)

# --- bacterial surface display (eCPX) --------------------------------------
BACTERIAL_SURFACE = Platform(
    key="bacteria_ecpx",
    name="Bacterial surface display (eCPX outer-membrane anchor)",
    name_zh="细菌表面展示（eCPX外膜锚定）",
    linkage="peptide on the outer membrane of the cell carrying its plasmid",
    compartment="E. coli outer membrane",
    chemistries=("disulfide", "linear"),
    noncanonical=False,
    readout="FACS against fluorescent target, with per-cell display normalisation",
    dna_input=1e12,
    practical_library=1e9,
    typical_ring_residues=(5, 20),
    steps=(
        Step("DNA", "Cassette ligation into the display vector",
             "盒式序列连入展示载体", "DNA",
             bulk=0.35,
             mechanism="eCPX terminal fusion cassette",
             consequence="bulk loss"),
        Step("XFORM", "Transformation", "转化", "host",
             bulk=0.05, ceiling=1e9,
             mechanism="the same electroporation ceiling as phage, usually "
                       "worked an order of magnitude below it",
             consequence="diversity ceiling"),
        Step("EXP", "Sec export and outer-membrane insertion",
             "Sec转运与外膜插入", "protein",
             bulk=0.55, sequence_term=sec_export,
             mechanism="the anchor must cross the periplasm and fold into the "
                       "outer membrane with the passenger attached",
             consequence="export bias, compounded by the anchor's own folding "
                         "requirement"),
        Step("FOLD", "Surface disulfide formation", "表面二硫键形成", "chemistry",
             sequence_term=lambda s: disulfide_closure(s)["efficiency"],
             mechanism="oxidising extracytoplasmic environment",
             consequence="requires a cysteine pair"),
        Step("PROT", "Surface proteolysis", "表面蛋白酶降解", "protein",
             sequence_term=proteolysis_survival,
             mechanism="OmpT is an outer-membrane protease acting directly on "
                       "the displayed passenger",
             consequence="dibasic motifs are depleted"),
        Step("TOX", "Growth-rate bias", "生长速率偏倚", "host",
             sequence_term=lambda s: 1.0 - host_toxicity(s),
             mechanism="membrane-active passengers burden the host",
             consequence="compounds over rounds"),
        Step("SORT", "FACS throughput", "流式分选通量", "readout",
             bulk=0.80, ceiling=1e8,
             mechanism="around 10^8 cells per sorting round is the practical "
                       "limit",
             consequence="quantitative per-cell readout bought at the cost of "
                         "throughput"),
    ),
)

# --- yeast display (the eukaryotic contrast, already a LitRaPID backend) ---
YEAST_DISPLAY = Platform(
    key="yeast_aga2",
    name="Yeast display (Aga2p fusion)",
    name_zh="酵母展示（Aga2p融合）",
    linkage="peptide on the cell wall, plasmid inside the same cell",
    compartment="yeast secretory pathway and cell wall",
    chemistries=("disulfide", "linear"),
    noncanonical=False,
    readout="two-colour FACS normalised for display level",
    dna_input=1e12,
    practical_library=1e7,
    typical_ring_residues=(6, 20),
    steps=(
        Step("DNA", "Gap-repair assembly", "同源重组组装", "DNA",
             bulk=0.40, mechanism="homologous recombination in vivo",
             consequence="bulk loss"),
        Step("XFORM", "Yeast transformation", "酵母转化", "host",
             bulk=0.01, ceiling=1e7,
             mechanism="~10^7 is the practical yeast library ceiling, two to "
                       "three orders below E. coli",
             consequence="the price of eukaryotic quality control is diversity"),
        Step("SEC", "ER translocation and quality control", "内质网转运与质量控制",
             "protein",
             bulk=0.45, sequence_term=sec_export,
             mechanism="chaperone and ERAD proofreading discards misfolded "
                       "passengers more strictly than the bacterial periplasm",
             consequence="cleaner display, smaller library"),
        Step("GLYC", "N-glycosylation of sequons", "糖基化修饰", "protein",
             sequence_term=glycosylation_survival,
             mechanism="N-X-S/T sequons are glycosylated in the secretory "
                       "pathway",
             consequence="an N-X-S/T ring is displayed as a different "
                         "chemical entity than the one designed"),
        Step("FOLD", "Oxidative folding in the ER", "内质网氧化折叠", "chemistry",
             sequence_term=lambda s: disulfide_closure(s)["efficiency"],
             mechanism="PDI-catalysed with chaperone proofreading",
             consequence="requires a cysteine pair"),
        Step("WALL", "Aga1p-Aga2p cell-wall anchoring", "Aga1p-Aga2p细胞壁锚定",
             "linkage",
             bulk=0.60, mechanism="disulfide-linked anchoring at 10^4-10^5 "
                                  "copies per cell",
             consequence="high avidity, which weakens affinity discrimination"),
        Step("SORT", "FACS throughput", "流式分选通量", "readout",
             bulk=0.80, ceiling=1e7, mechanism="sorting throughput limit",
             consequence="ceiling"),
    ),
)

PLATFORMS: dict[str, Platform] = {
    p.key: p for p in (MRNA_DISPLAY, PHAGE_DISPLAY, BACTERIAL_SICLOPPS,
                       BACTERIAL_SURFACE, YEAST_DISPLAY)
}


# ---------------------------------------------------------------------------
# Per-sequence profile
# ---------------------------------------------------------------------------
def biosynthesis_profile(sequence: str, platform: str | Platform) -> dict[str, Any]:
    """Walk one sequence through one platform's route.

    Returns the per-step survival probability, the overall probability that a
    given template yields a correctly displayed macrocycle, and the step that
    limits it. The limiting step is the actionable output: it names the
    experiment to change, which an aggregate score cannot.
    """
    plat = platform if isinstance(platform, Platform) else PLATFORMS[platform]
    seq = sequence.strip()
    if not seq:
        raise ValueError("sequence must not be empty")
    unknown = {r for r in seq if r not in AA and r not in NONCANONICAL and r != "*"}
    if unknown:
        raise ValueError(f"unknown residue code(s): {''.join(sorted(unknown))}")
    if not plat.noncanonical and any(r in NONCANONICAL for r in seq):
        return {
            "platform": plat.key, "sequence": seq, "display_probability": 0.0,
            "limiting_step": "genetic_code",
            "limiting_reason": "non-canonical residues cannot be encoded by a "
                               "cell-based platform; this molecule requires "
                               "flexizyme-charged tRNAs",
            "steps": [], "chemistry": None,
        }

    current = seq
    steps: list[dict[str, Any]] = []
    total = 1.0
    for step in plat.steps:
        if step.rewrites is not None:
            current = step.rewrites(current)
        term = step.sequence_term(current) if step.sequence_term else 1.0
        survival = step.bulk * term
        total *= survival
        steps.append({
            "token": step.token, "name": step.name, "name_zh": step.name_zh,
            "layer": step.layer, "category": step.category,
            "bulk": round(step.bulk, 6),
            "sequence_term": round(term, 6),
            "survival": round(survival, 8),
            "cumulative": round(total, 10),
            "mechanism": step.mechanism, "consequence": step.consequence,
        })

    informative = [s for s in steps if s["survival"] < 0.999]
    limiting = min(informative, key=lambda s: s["survival"]) if informative else steps[0]

    chemistry: dict[str, Any] | None = None
    if "thioether" in plat.chemistries:
        chemistry = {"route": "thioether", **thioether_closure(current)}
    elif plat.key == "bacteria_siclopps":
        chemistry = {"route": "head_to_tail", **intein_splicing(current, "Npu")}
    elif "disulfide" in plat.chemistries:
        chemistry = {"route": "disulfide", **disulfide_closure(current)}

    return {
        "platform": plat.key,
        "platform_name": plat.name,
        "sequence": seq,
        "expressed_sequence": current,
        "display_probability": total,
        "limiting_step": limiting["token"],
        "limiting_reason": limiting["consequence"],
        "steps": steps,
        "chemistry": chemistry,
        "descriptors": {
            "length": len(seq),
            "net_charge": net_charge(seq),
            "hydropathy": round(hydropathy(seq), 4),
            "cysteines": seq.count("C"),
            "n_methyl": n_methyl_count(seq),
            "dibasic_sites": sum(1 for i in range(len(seq) - 1)
                                 if seq[i] in BASIC and seq[i + 1] in BASIC),
        },
    }


# ---------------------------------------------------------------------------
# Library-scale funnel
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class LibrarySpec:
    """The randomised design handed to the platform."""

    length: int = 12
    codon_design: str = "NNK"
    fixed: dict[int, str] = field(default_factory=dict)
    input_molecules: Optional[float] = None
    sample_size: int = 4000
    seed: int = 2025

    def validate(self) -> None:
        if not 3 <= self.length <= 40:
            raise ValueError("length must be between 3 and 40")
        if self.codon_design not in STOP_FREQUENCY:
            raise ValueError(f"unsupported codon design: {self.codon_design}")
        if not 100 <= self.sample_size <= 200_000:
            raise ValueError("sample_size must be between 100 and 200,000")
        if self.input_molecules is not None and self.input_molecules <= 0:
            raise ValueError("input_molecules must be positive")

    @property
    def nominal_diversity(self) -> float:
        return 20.0 ** (self.length - len(self.fixed))


def _codon_residue_weights(codon_design: str) -> dict[str, float]:
    """Amino-acid frequencies implied by a degenerate codon, via the code.

    A display library is synthesised as DNA, so the residue distribution is the
    image of a nucleotide distribution under the genetic code and is never
    uniform: with NNK, Leu/Arg/Ser get three codons each while Met/Trp get one.
    Sampling residues uniformly models a library that exists in no tube.
    """
    if codon_design == "trimer":
        return {a: 1.0 / 20 for a in AA}
    codons = _expand_codons(codon_design)
    weights: dict[str, float] = {}
    for codon in codons:
        residue = CODON_TABLE[codon]
        weights[residue] = weights.get(residue, 0.0) + 1.0
    return {k: v / len(codons) for k, v in weights.items()}


def sample_library(spec: LibrarySpec, platform: Platform,
                   rng: random.Random) -> list[str]:
    """Draw representative members, respecting the platform's fixed chemistry.

    SICLOPPS fixes a cysteine at the extein's first position; RaPID libraries
    conventionally fix a cysteine at the C-terminal end of the randomised
    region so the chloroacetamide has a partner. Those fixed residues are part
    of the platform, not of the design, and omitting them makes the cyclisation
    statistics meaningless.
    """
    weights = _codon_residue_weights(spec.codon_design)
    residues = list(weights)
    probabilities = [weights[r] for r in residues]
    fixed = dict(spec.fixed)
    if platform.key == "bacteria_siclopps":
        fixed.setdefault(0, "C")
    elif platform.key == "mrna_rapid":
        fixed.setdefault(spec.length - 1, "C")

    members: list[str] = []
    for _ in range(spec.sample_size):
        draw = rng.choices(residues, weights=probabilities, k=spec.length)
        for position, residue in fixed.items():
            if 0 <= position < spec.length:
                draw[position] = residue
        sequence = "".join(draw)
        if "*" in sequence:
            # a stop truncates the open reading frame at that codon
            sequence = sequence[:sequence.index("*")] + "*"
        members.append(sequence)
    return members


def _amber_fate(platform: Platform, sequence: str) -> tuple[str, float]:
    """What a platform does with the amber codon an NNK design inevitably puts in.

    This is where the three routes genuinely differ. In a PURE system lacking
    RF1 the codon is reassignable, so the member survives as a full-length
    product. In a supE host it is read as glutamine most of the time, so the
    member survives with its sequence silently rewritten. Everywhere else it
    terminates and the member is simply gone.
    """
    if "*" not in sequence:
        return sequence, 1.0
    if platform.noncanonical:                       # PURE, no RF1
        return sequence.replace("*", "B"), 0.55     # reassigned to an ncAA
    if platform.key == "phage_m13":                 # supE suppression
        return sequence.replace("*", "Q"), 0.60
    return sequence.replace("*", ""), 0.0           # terminates


def simulate_biosynthesis(platform: str | Platform,
                          spec: LibrarySpec | None = None) -> dict[str, Any]:
    """Propagate a library through one platform's biosynthetic route.

    Molecules and distinct variants are propagated separately. For a variant
    present in n copies and surviving a step with probability p, the chance
    that at least one copy survives is 1 - exp(-n*p); summing that over the
    sampled sequence distribution gives the distinct count, while summing n*p
    gives the molecule count. The gap between the two curves is the difference
    between losing material and losing library.
    """
    plat = platform if isinstance(platform, Platform) else PLATFORMS[platform]
    spec = spec or LibrarySpec()
    spec.validate()
    rng = random.Random(spec.seed)

    members = sample_library(spec, plat, rng)
    fates = [_amber_fate(plat, m) for m in members]

    dna_input = float(spec.input_molecules or plat.dna_input)
    molecules = dna_input
    distinct = min(spec.nominal_diversity, molecules)

    # The amber codon is resolved before the route begins: it is a property of
    # the template and the translation system, not of a later physical step.
    amber_survival = sum(f[1] for f in fates) / len(fates)
    sequences = [f[0] for f in fates]
    weights = [f[1] for f in fates]

    rows: list[dict[str, Any]] = []
    rows.append({
        "token": "LIB", "name": "Library design as synthesised",
        "name_zh": "设计文库（合成态）", "layer": "DNA", "category": "bulk",
        "molecules": molecules, "distinct": distinct,
        "step_survival": 1.0, "distinct_retained": 1.0,
        "mechanism": f"{spec.codon_design} randomised at {spec.length} positions",
        "consequence": f"nominal diversity {spec.nominal_diversity:.3g}",
    })

    copies = molecules / distinct if distinct else 0.0
    distinct_after = sum(1.0 - math.exp(-copies * w) for w in weights) / len(weights)
    molecules *= amber_survival
    distinct *= distinct_after
    rows.append({
        "token": "AMBER", "name": "Fate of the amber codon",
        "name_zh": "琥珀密码子的去向", "layer": "translation", "category": "sequence",
        "molecules": molecules, "distinct": distinct,
        "step_survival": amber_survival, "distinct_retained": distinct_after,
        "mechanism": "NNK places one stop codon per 32; "
                     + ("RF1 is absent, so it is reassignable"
                        if plat.noncanonical else
                        "supE reads it as Gln" if plat.key == "phage_m13"
                        else "it terminates translation"),
        "consequence": f"stop-free fraction of the design is "
                       f"{stop_free_fraction(spec.codon_design, spec.length):.3f}",
    })

    cumulative = [w for w in weights]
    for step in plat.steps:
        if step.rewrites is not None:
            sequences = [step.rewrites(s) for s in sequences]
        terms = [step.bulk * (step.sequence_term(s) if step.sequence_term else 1.0)
                 for s in sequences]
        mean_survival = sum(terms) / len(terms)

        copies = molecules / distinct if distinct > 0 else 0.0
        retained = sum(1.0 - math.exp(-copies * t) for t in terms) / len(terms)

        molecules *= mean_survival
        distinct *= retained
        cumulative = [c * t for c, t in zip(cumulative, terms)]

        capped = False
        if step.ceiling is not None and distinct > step.ceiling:
            molecules *= step.ceiling / distinct
            distinct = step.ceiling
            capped = True

        rows.append({
            "token": step.token, "name": step.name, "name_zh": step.name_zh,
            "layer": step.layer, "category": step.category,
            "molecules": molecules, "distinct": distinct,
            "step_survival": mean_survival, "distinct_retained": retained,
            "ceiling_applied": capped,
            "mechanism": step.mechanism, "consequence": step.consequence,
        })

    naive = _residue_frequencies(members)
    displayed = _residue_frequencies(sequences, cumulative)
    bias = {
        residue: round(math.log2((displayed.get(residue, 0.0) + 1e-6)
                                 / (naive.get(residue, 0.0) + 1e-6)), 4)
        for residue in AA
    }

    informative = [r for r in rows[1:] if r["step_survival"] < 0.999]
    diversity_losses = [r for r in rows[1:] if r["distinct_retained"] < 0.999
                        or r.get("ceiling_applied")]

    return {
        "platform": plat.key,
        "platform_name": plat.name,
        "platform_name_zh": plat.name_zh,
        "linkage": plat.linkage,
        "compartment": plat.compartment,
        "chemistries": list(plat.chemistries),
        "noncanonical_capable": plat.noncanonical,
        "library": asdict(spec),
        "funnel": rows,
        "final_molecules": rows[-1]["molecules"],
        "final_distinct": rows[-1]["distinct"],
        "dna_input_molecules": dna_input,
        "overall_yield": rows[-1]["molecules"] / dna_input,
        "diversity_retained": rows[-1]["distinct"] / spec.nominal_diversity,
        "yield_limiting_step": min(informative,
                                   key=lambda r: r["step_survival"])["token"]
                               if informative else None,
        "diversity_limiting_step": min(diversity_losses,
                                       key=lambda r: r["distinct"])["token"]
                                   if diversity_losses else None,
        "residue_bias_log2": bias,
        "mean_display_probability": sum(cumulative) / len(cumulative),
        "limitations": "Mechanistically-motivated parameterisation, not a "
                       "fitted model. Use for cross-platform and "
                       "cross-sequence comparison; recalibrate per-step "
                       "constants on measured data before quoting yields.",
    }


def _residue_frequencies(sequences: list[str],
                         weights: list[float] | None = None) -> dict[str, float]:
    weights = weights or [1.0] * len(sequences)
    counts: dict[str, float] = {}
    total = 0.0
    for sequence, weight in zip(sequences, weights):
        if weight <= 0:
            continue
        for residue in sequence:
            key = residue_properties(residue)[0]
            counts[key] = counts.get(key, 0.0) + weight
            total += weight
    return {k: v / total for k, v in counts.items()} if total else {}


def compare_platforms(spec: LibrarySpec | None = None,
                      platforms: list[str] | None = None) -> dict[str, Any]:
    """Run the same library design through every platform and tabulate the result."""
    keys = platforms or list(PLATFORMS)
    results = {k: simulate_biosynthesis(k, spec) for k in keys}
    table = [{
        "platform": k,
        "name": r["platform_name"],
        "chemistries": ",".join(r["chemistries"]),
        "noncanonical": r["noncanonical_capable"],
        "final_distinct": r["final_distinct"],
        "diversity_retained": r["diversity_retained"],
        "overall_yield": r["overall_yield"],
        "yield_limiting_step": r["yield_limiting_step"],
        "diversity_limiting_step": r["diversity_limiting_step"],
    } for k, r in results.items()]
    table.sort(key=lambda row: row["final_distinct"], reverse=True)
    return {"library": asdict(spec or LibrarySpec()), "comparison": table,
            "detail": results}


def display_probability_table(sequences: list[str],
                              platform: str | Platform) -> dict[str, float]:
    """Per-sequence display probabilities, for injection into a selection model.

    `mrna_display.simulate_mrna_display` multiplies every candidate by one
    shared `fusion_probability`. Passing this table instead makes the
    biosynthetic prior sequence-resolved, which is what lets a selection
    simulation attribute enrichment to binding rather than to manufacturability.
    """
    out: dict[str, float] = {}
    for sequence in sequences:
        try:
            out[sequence] = biosynthesis_profile(sequence, platform)["display_probability"]
        except ValueError:
            out[sequence] = 0.0
    return out


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="litrapid.biosynthesis",
        description="Mechanistic biosynthesis simulation of displayed macrocycles")
    sub = parser.add_subparsers(dest="command", required=True)

    profile = sub.add_parser("profile", help="walk one sequence through one platform")
    profile.add_argument("--sequence", required=True)
    profile.add_argument("--platform", required=True, choices=sorted(PLATFORMS))

    funnel = sub.add_parser("funnel", help="library-scale funnel for one platform")
    funnel.add_argument("--platform", required=True, choices=sorted(PLATFORMS))
    funnel.add_argument("--length", type=int, default=12)
    funnel.add_argument("--codon-design", default="NNK", choices=sorted(STOP_FREQUENCY))
    funnel.add_argument("--sample-size", type=int, default=4000)
    funnel.add_argument("--seed", type=int, default=2025)

    compare = sub.add_parser("compare", help="all platforms on the same design")
    compare.add_argument("--length", type=int, default=12)
    compare.add_argument("--codon-design", default="NNK", choices=sorted(STOP_FREQUENCY))
    compare.add_argument("--sample-size", type=int, default=4000)
    compare.add_argument("--seed", type=int, default=2025)
    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    if args.command == "profile":
        payload: Any = biosynthesis_profile(args.sequence, args.platform)
    else:
        spec = LibrarySpec(length=args.length, codon_design=args.codon_design,
                           sample_size=args.sample_size, seed=args.seed)
        if args.command == "funnel":
            payload = simulate_biosynthesis(args.platform, spec)
        else:
            payload = compare_platforms(spec)
            payload.pop("detail")
    print(json.dumps(payload, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
