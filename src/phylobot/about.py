"""The About page's citations: every method PhyloBot runs, and the substitution models.

Mirrors the README's "Included methods" section; a change to one belongs in the other.
Every DOI here was checked against Crossref (2026-10-01). Versions are those pinned in
pipeline/nextflow.config.
"""

import csv
from dataclasses import dataclass
from pathlib import Path

from phylobot import catalog


@dataclass(frozen=True)
class Ref:
    authors: str
    where: str  # journal, volume, pages, year
    doi: str = ""
    title: str = ""
    url: str = ""  # for references without a DOI

    @property
    def link(self) -> str:
        return f"https://doi.org/{self.doi}" if self.doi else self.url


@dataclass(frozen=True)
class Method:
    name: str
    note: str
    refs: tuple[Ref, ...] = ()


PHYLOBOT_2016 = Ref(
    "Hanson-Smith V, Johnson A.",
    "PLoS Computational Biology 12(7): e1004976 (2016).",
    "10.1371/journal.pcbi.1004976",
    "PhyloBot: A Web Portal for Automated Phylogenetics, Ancestral Sequence "
    "Reconstruction, and Exploration of Mutational Trajectories.",
)

_PRANK = (
    Ref(
        "Löytynoja A, Goldman N.",
        "PNAS 102(30): 10557–10562 (2005).",
        "10.1073/pnas.0409137102",
        "An algorithm for progressive multiple alignment of sequences with insertions.",
    ),
    Ref(
        "Löytynoja A, Goldman N.",
        "Science 320(5883): 1632–1635 (2008).",
        "10.1126/science.1158395",
        "Phylogeny-aware gap placement prevents errors in sequence alignment and "
        "evolutionary analysis.",
    ),
)
_VANANTWERP = Ref(
    "VanAntwerp J, Finneran P, Dolgikh B, Woldring D.",
    "Methods in Molecular Biology 2491: 75–86 (2022).",
    "10.1007/978-1-0716-2285-8_4",
    "Ancestral Sequence Reconstruction and Alternate Amino Acid States Guide Protein "
    "Library Design for Directed Evolution.",
)

_POSADA = Ref(
    "Posada D, Buckley TR.",
    "Systematic Biology 53(5): 793–808 (2004).",
    "10.1080/10635150490522304",
    "Model Selection and Model Averaging in Phylogenetics: Advantages of Akaike "
    "Information Criterion and Bayesian Approaches over Likelihood Ratio Tests.",
)

_GONZALEZ = Ref(
    "Gonzalez TF.",
    "Theoretical Computer Science 38: 293–306 (1985).",
    "10.1016/0304-3975(85)90224-5",
    "Clustering to minimize the maximum intercluster distance.",
)

_EICK = Ref(
    "Eick GN, Bridgham JT, Anderson DP, Harms MJ, Thornton JW.",
    "Molecular Biology and Evolution 34(2): 247–261 (2017).",
    "10.1093/molbev/msw223",
    "Robustness of Reconstructed Ancestral Protein Functions to Statistical Uncertainty.",
)

_IQTREE = Ref(
    "Wong TKF, Ly-Trong N, Ren H, et al.",
    "Molecular Biology and Evolution 43(5): msag117 (2026).",
    "10.1093/molbev/msag117",
    "IQ-TREE 3: phylogenomic inference software using complex evolutionary models.",
)

# Sequence properties measured over an ancestor's posterior. The scales
# are the published ones; the summaries over them (windows, thresholds, the charge model)
# are PhyloBot's own and are described, not cited, so nobody reads a citation as endorsing
# our particular choice of window or cutoff.
_EISENBERG_MOMENT = Ref(
    "Eisenberg D, Weiss RM, Terwilliger TC.",
    "Nature 299: 371–374 (1982).",
    "10.1038/299371a0",
    "The helical hydrophobic moment: a measure of the amphiphilicity of a helix.",
)
_EISENBERG_SCALE = Ref(
    "Eisenberg D, Schwarz E, Komaromy M, Wall R.",
    "Journal of Molecular Biology 179(1): 125–142 (1984).",
    "10.1016/0022-2836(84)90309-7",
    "Analysis of membrane and surface protein sequences with the hydrophobic moment plot.",
)
_KYTE_DOOLITTLE = Ref(
    "Kyte J, Doolittle RF.",
    "Journal of Molecular Biology 157(1): 105–132 (1982).",
    "10.1016/0022-2836(82)90515-0",
    "A simple method for displaying the hydropathic character of a protein.",
)

# The neural predictors, run in the phylobot-ml image. DOIs checked
# against Crossref, 2026-10-21; metapredict V3 is a preprint.
_METAPREDICT = Ref(
    "Emenecker RJ, Griffith D, Holehouse AS.",
    "Biophysical Journal 120(20): 4312–4319 (2021).",
    "10.1016/j.bpj.2021.08.039",
    "Metapredict: a fast, accurate, and easy-to-use predictor of consensus disorder and structure.",
)
_METAPREDICT_V3 = Ref(
    "Lotthammer JM, Hernández-García J, Griffith D, Weijers D, Holehouse AS, Emenecker RJ.",
    "bioRxiv 2024.11.05.622168 (2024), preprint.",
    "10.1101/2024.11.05.622168",
    "Metapredict enables accurate disorder prediction across the Tree of Life.",
)
_ESM2 = Ref(
    "Lin Z, Akin H, Rao R, Hie B, Zhu Z, Lu W, et al.",
    "Science 379(6637): 1123–1130 (2023).",
    "10.1126/science.ade2574",
    "Evolutionary-scale prediction of atomic-level protein structure with a language model.",
)


METHODS: tuple[tuple[str, tuple[Method, ...]], ...] = (
    ("Multiple sequence alignment", (
        Method("MAFFT", "7.525", (Ref(
            "Katoh K, Standley DM.",
            "Molecular Biology and Evolution 30(4): 772–780 (2013).",
            "10.1093/molbev/mst010",
            "MAFFT Multiple Sequence Alignment Software Version 7: Improvements in "
            "Performance and Usability.",
        ),)),
        Method("MUSCLE", "5.3", (Ref(
            "Edgar RC.",
            "Nature Communications 13: 6968 (2022).",
            "10.1038/s41467-022-34630-w",
            "Muscle5: High-accuracy alignment ensembles enable unbiased assessments of "
            "sequence homology and phylogeny.",
        ),)),
        Method("Clustal Omega", "1.2.4", (Ref(
            "Sievers F, Wilm A, Dineen D, et al.",
            "Molecular Systems Biology 7: 539 (2011).",
            "10.1038/msb.2011.75",
            "Fast, scalable generation of high-quality protein multiple sequence "
            "alignments using Clustal Omega.",
        ),)),
        Method("PRANK", "v.250331; phylogeny-aware, keeps independent insertions apart",
               _PRANK),
    )),
    ("Phylogeny, model selection and branch support", (
        Method("IQ-TREE", "3.1.3; tree search, model fitting, ancestral reconstruction",
               (_IQTREE,)),
        Method("RAxML-NG", "2.0.3; a second tree search, whose trees IQ-TREE then fits", (
            Ref(
                "Kozlov AM, Darriba D, Flouri T, Morel B, Stamatakis A.",
                "Bioinformatics 35(21): 4453–4455 (2019).",
                "10.1093/bioinformatics/btz305",
                "RAxML-NG: a fast, scalable and user-friendly tool for maximum "
                "likelihood phylogenetic inference.",
            ),
            Ref(
                "Kozlov O, Togkousidis A, Stelz C, et al.",
                "bioRxiv preprint (2026).",
                "10.64898/2026.09.09.750097",
                "RAxML-NG 2: Automatic model selection, novel tree search heuristics, "
                "and fast branch support metrics.",
            ),
        )),
        Method("GeneRax", "2.1.3; reconciles each tree with the species tree, and "
               "searches a reconciled tree", (Ref(
            "Morel B, Kozlov AM, Stamatakis A, Szöllősi GJ.",
            "Molecular Biology and Evolution 37(9): 2763–2774 (2020).",
            "10.1093/molbev/msaa141",
            "GeneRax: a tool for species-tree-aware maximum likelihood-based gene family "
            "tree inference under gene duplication, transfer, and loss.",
        ),)),
        Method("Open Tree of Life", "the species tree for reconciliation (synthetic tree and "
               "taxonomy; the release is recorded with each job)", (
            Ref(
                "Hinchliff CE, Smith SA, Allman JF, et al.",
                "Proceedings of the National Academy of Sciences 112(41): 12764–12769 (2015).",
                "10.1073/pnas.1423041112",
                "Synthesis of phylogeny and taxonomy into a comprehensive tree of life.",
            ),
            Ref(
                "Redelings BD, Holder MT.",
                "PeerJ 5: e3058 (2017).",
                "10.7717/peerj.3058",
                "A supertree pipeline for summarizing phylogenetic and taxonomic information "
                "for millions of species.",
            ),
        )),
        Method("UniProt", "species read from UniProt entry names (its controlled vocabulary "
               "of species, CC BY 4.0, bundled) and from UniProt accessions, with consent", (Ref(
            "The UniProt Consortium.",
            "Nucleic Acids Research 53(D1): D609–D617 (2025).",
            "10.1093/nar/gkae1010",
            "UniProt: the Universal Protein Knowledgebase in 2025.",
        ),)),
        Method("AU test", "in IQ-TREE; is each tree source's topology rejected by the "
               "sequences?", (Ref(
            "Shimodaira H.",
            "Systematic Biology 51(3): 492–508 (2002).",
            "10.1080/10635150290069913",
            "An approximately unbiased test of phylogenetic tree selection.",
        ),)),
        Method("ModelFinder", "in IQ-TREE", (Ref(
            "Kalyaanamoorthy S, Minh BQ, Wong TKF, von Haeseler A, Jermiin LS.",
            "Nature Methods 14: 587–589 (2017).",
            "10.1038/nmeth.4285",
            "ModelFinder: fast model selection for accurate phylogenetic estimates.",
        ),)),
        Method("Model comparison", "AIC, AICc and BIC, with weights", (
            Ref("Akaike H.", "IEEE Transactions on Automatic Control 19(6): 716–723 (1974).",
                "10.1109/TAC.1974.1100705", "A new look at the statistical model identification."),
            Ref("Hurvich CM, Tsai CL.", "Biometrika 76(2): 297–307 (1989).",
                "10.1093/biomet/76.2.297",
                "Regression and time series model selection in small samples."),
            Ref("Schwarz G.", "Annals of Statistics 6(2): 461–464 (1978).",
                "10.1214/aos/1176344136", "Estimating the Dimension of a Model."),
        )),
        Method("UFBoot", "ultrafast bootstrap", (Ref(
            "Hoang DT, Chernomor O, von Haeseler A, Minh BQ, Vinh LS.",
            "Molecular Biology and Evolution 35(2): 518–522 (2018).",
            "10.1093/molbev/msx281",
            "UFBoot2: Improving the Ultrafast Bootstrap Approximation.",
        ),)),
        Method("SH-aLRT", "SH-like approximate likelihood ratio test", (Ref(
            "Guindon S, Dufayard JF, Lefort V, Anisimova M, Hordijk W, Gascuel O.",
            "Systematic Biology 59(3): 307–321 (2010).",
            "10.1093/sysbio/syq010",
            "New Algorithms and Methods to Estimate Maximum-Likelihood Phylogenies: "
            "Assessing the Performance of PhyML 3.0.",
        ),)),
    )),
    ("Ancestral sequence reconstruction", (
        Method("Empirical Bayes (marginal) reconstruction", "IQ-TREE --ancestral", (Ref(
            "Yang Z, Kumar S, Nei M.",
            "Genetics 141(4): 1641–1650 (1995).",
            "10.1093/genetics/141.4.1641",
            "A new method of inference of ancestral nucleotide and amino acid sequences.",
        ),)),
    )),
    ("Gaps in ancestral sequences", (
        Method("Parsimony (shared presence)",
               "PhyloBot's own rule: a residue only where it is present on both sides "
               "of a node"),
        Method("Parsimony (Fitch)", "Hartigan's rule for multifurcations", (
            Ref("Fitch WM.", "Systematic Zoology 20(4): 406–416 (1971).", "10.2307/2412116",
                "Toward Defining the Course of Evolution: Minimum Change for a Specific "
                "Tree Topology."),
            Ref("Hartigan JA.", "Biometrics 29(1): 53–65 (1973).", "10.2307/2529676",
                "Minimum Mutation Fits to a Given Tree."),
        )),
        Method("Maximum likelihood",
               "the gap/residue pattern reconstructed by IQ-TREE under a two-state model "
               "on the fixed tree", (_IQTREE,)),
        Method("PRANK", "its own indel-aware ancestors (-showanc)", _PRANK),
        Method("GRASP", "2.0.0; joint reconstruction, bi-directional edge parsimony for "
               "indels", (Ref(
                   "Foley G, Mora A, Ross CM, et al.",
                   "PLoS Computational Biology 18(10): e1010633 (2022).",
                   "10.1371/journal.pcbi.1010633",
                   "Engineering indel and substitution variants of diverse and ancient "
                   "enzymes using Graphical Representation of Ancestral Sequence "
                   "Predictions (GRASP).",
               ),)),
    )),
    # The ancestral library. AltAll and the single-site alternates come from the
    # literature; the three fill strategies are PhyloBot's own, so each says what it is
    # and what inspired it rather than claiming a method paper. They were measured on a
    # simulated family (2026-10-13).
    ("Ancestral libraries", (
        Method("AltAll", "one sequence carrying the alternative state at every ambiguous "
               "site: the ASR literature's robustness control, and a library's second "
               "well", (_EICK,)),
        Method("Posterior samples", "draws from IQ-TREE's per-site marginal posteriors "
               "(no joint distribution is available), which Eick et al. use alongside "
               "AltAll to explore plausible sequences", (_EICK, _IQTREE)),
        Method("Single-site alternates", "one well per ambiguous site, or per pair: "
               "site-wise diversification around the ancestor", (_VANANTWERP,)),
        Method("Most supported", "PhyloBot's own: the sequences that the most distinct "
               "choices of aligner, model and gap placement agree on, ranked by their "
               "summed BIC weight. The inspiration is model averaging, where a result "
               "is weighted by how well each model fits rather than taken from one",
               (_POSADA,)),
        Method("Spanning", "PhyloBot's own: greedy max–min (farthest-first) selection, so "
               "the wells cover the disagreement between reconstructions instead of "
               "holding near-identical variants. The algorithm is the standard k-centre "
               "greedy traversal; its use here is not a published ASR method, and on "
               "PhyloBot's own simulated test it did not improve accuracy over Most "
               "supported", (_GONZALEZ,)),
        Method("Per ancestor", "PhyloBot's own: equal shares of the wells per ancestor, "
               "round robin, for a library over a lineage. An allocation rule, with no "
               "published precedent and no statistical claim"),
    )),
    ("Sequence properties", (
        Method(
            "Physicochemical descriptors",
            "measured over an ancestral posterior: net charge, the hydrophobic "
            "moment on the Eisenberg scale and hydropathy runs on Kyte–Doolittle, "
            "summarised by PhyloBot's own windows and thresholds",
            (_EISENBERG_MOMENT, _EISENBERG_SCALE, _KYTE_DOOLITTLE),
        ),
        Method("metapredict", "3.0.2 (V3 network): the disordered fraction of an "
               "ancestor, optional", (_METAPREDICT, _METAPREDICT_V3)),
        Method("ESM-2", "35M parameters (esm2_t12_35M_UR50D): how natural an ancestor "
               "looks to a protein language model, optional", (_ESM2,)),
    )),
    ("Homolog search", (
        Method("BLAST", "blastp", (Ref(
            "Camacho C, Coulouris G, Avagyan V, et al.",
            "BMC Bioinformatics 10: 421 (2009).",
            "10.1186/1471-2105-10-421",
            "BLAST+: architecture and applications.",
        ),)),
        Method("EMBL-EBI Job Dispatcher", "remote searches at EBI", (Ref(
            "Madeira F, Madhusoodanan N, Lee J, et al.",
            "Nucleic Acids Research 52(W1): W521–W525 (2024).",
            "10.1093/nar/gkae241",
            "The EMBL-EBI Job Dispatcher sequence analysis tools framework in 2024.",
        ),)),
        Method("UniRef90", "the database searched at EBI", (Ref(
            "Suzek BE, Wang Y, Huang H, McGarvey PB, Wu CH, UniProt Consortium.",
            "Bioinformatics 31(6): 926–932 (2015).",
            "10.1093/bioinformatics/btu739",
            "UniRef clusters: a comprehensive and scalable alternative for improving "
            "sequence similarity searches.",
        ),)),
        Method("NCBI BLAST", "ClusteredNR; local mode only", (Ref(
            "National Center for Biotechnology Information.", "BLAST web service.",
            url="https://blast.ncbi.nlm.nih.gov/",
        ),)),
    )),
    ("Workflow and testing", (
        Method("Nextflow", "runs the pipeline", (Ref(
            "Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C.",
            "Nature Biotechnology 35: 316–319 (2017).",
            "10.1038/nbt.3820",
            "Nextflow enables reproducible computational workflows.",
        ),)),
        Method("BioContainers", "the pinned tool images", (Ref(
            "da Veiga Leprevost F, Grüning BA, Alves Aflitos S, et al.",
            "Bioinformatics 33(16): 2580–2582 (2017).",
            "10.1093/bioinformatics/btx192",
            "BioContainers: an open-source and community-driven framework for software "
            "standardization.",
        ),)),
        Method("DendroPy", "tree handling in PhyloBot's own steps", (Ref(
            "Moreno MA, Holder MT, Sukumaran J.",
            "Journal of Open Source Software 9(101): 6943 (2024).",
            "10.21105/joss.06943",
            "DendroPy 5: a mature Python library for phylogenetic computing.",
        ),)),
        Method("AliSim", "the simulated test datasets with known ancestors", (Ref(
            "Ly-Trong N, Naser-Khdour S, Lanfear R, Minh BQ.",
            "Molecular Biology and Evolution 39(5): msac092 (2022).",
            "10.1093/molbev/msac092",
            "AliSim: A Fast and Versatile Phylogenetic Sequence Simulator for the "
            "Genomic Era.",
        ),)),
    )),
)  # fmt: skip


#: One entry per `phylobot.properties.REGISTRY` key: how the number is computed, in
#: enough detail to reproduce it, and what it is *not*. A test checks every registered
#: property has one.
PROPERTY_METHODS: dict[str, Method] = {
    "charge": Method(
        "Net charge (pH 7)",
        "Answers: is this ancestor overall positively or negatively charged, and by how "
        "much? Charge shapes who a protein binds — nucleic acids, membranes, other "
        "proteins — and whether it stays soluble. "
        "Computed as a sum of fixed per-residue charges over the ancestor's residues, "
        "gaps excluded: Lys and Arg +1, Asp and Glu −1, His +0.1 (its pKa is near 6, so "
        "it carries little charge at pH 7). This is a counting model, not a pKa "
        "calculation: it ignores the structure, local environment and terminal groups, "
        "so treat it as a composition descriptor rather than a predicted titration "
        "state.",
    ),
    "moment": Method(
        "Hydrophobic moment",
        "Answers: does this ancestor have stretches that are oily on one face and "
        "water-liking on the other? That lopsided pattern is the signature of a helix "
        "that sits at a surface — a membrane, an interface, a binding groove. "
        "Computed as the mean helical hydrophobic moment over sliding 11-residue "
        "windows, each residue placed at 100° of helical phase, on the Eisenberg "
        "consensus scale; reported per residue. It measures how far hydrophobicity "
        "follows the helical period, so a uniformly hydrophobic sequence scores near "
        "zero. The window length and the averaging over all windows are PhyloBot's "
        "choices.",
        (_EISENBERG_MOMENT, _EISENBERG_SCALE),
    ),
    "hydropathy_regions": Method(
        "Hydrophobic regions",
        "Answers: how many separate oily patches does this ancestor have? Such patches "
        "are what get buried in a folded core, cross a membrane, or — when exposed — "
        "stick to things they should not. "
        "Computed as the number of distinct runs of sliding 7-residue windows whose "
        "mean Kyte–Doolittle hydropathy exceeds 1.5. Deliberately non-smooth — one "
        "substitution can tip a window across the threshold and split or merge a run — "
        "because it stands in for a real aggregation predictor and tests whether the "
        "machinery survives a discontinuous property. The window and threshold are "
        "PhyloBot's; it is not a validated aggregation or transmembrane predictor.",
        (_KYTE_DOOLITTLE,),
    ),
    "disorder": Method(
        "Disordered fraction",
        "Answers: how much of this ancestor had no fixed shape? Disordered regions stay "
        "floppy instead of folding, and they are where much signalling, regulation and "
        "partner-switching happens — so a shift here is a shift in what kind of protein "
        "this was. "
        "Computed as the share of the ancestor's residues (gaps excluded) inside "
        "intrinsically disordered regions called by metapredict 3.0.2 with its V3 "
        "network and its own defaults (regions of at least 12 residues, gaps of up to "
        "10 closed, folded domains of at least 50). A sequence-only predictor trained "
        "on consensus disorder of natural proteins: on a reconstructed ancestor it says "
        "what the sequence resembles, not what the protein did, and a draw that is a "
        "worse protein than any real one can look more disordered for that reason "
        "alone.",
        (_METAPREDICT, _METAPREDICT_V3),
    ),
    "esm_likelihood": Method(
        "ESM-2 log-likelihood",
        "Answers: does this ancestor look like a real protein? A model that has read "
        "millions of natural sequences scores how unsurprising each residue is, which "
        "is a way of asking whether a reconstruction is plausible or an artefact. "
        "Higher is more natural. "
        "Computed as the mean, over the ancestor's residues (gaps excluded), of the "
        "natural-log probability ESM-2 35M assigns to each residue in one forward pass "
        "with nothing masked; proteins longer than its 1022-residue context are scored "
        "in overlapping windows. One unmasked pass stands in for the "
        "pseudo-log-likelihood (one masked pass per residue), which costs L times more; "
        "on Ime2 it ranked real sequences above 10%-mutated and shuffled ones, though "
        "not mutants against shuffles consistently. It is a plausibility score, not a "
        "fitness, stability or activity prediction, and its values compare only with "
        "each other: an ancestor's draws against its extant relatives.",
        (_ESM2,),
    ),
}


@dataclass(frozen=True)
class ModelRow:
    names: tuple[str, ...]  # matrix names, or variant suffixes
    use: str
    refs: tuple[Ref, ...] = ()


MATRICES: tuple[ModelRow, ...] = (
    ModelRow(("LG",), "general", (Ref(
        "Le SQ, Gascuel O.", "Mol Biol Evol 25(7): 1307–1320 (2008).",
        "10.1093/molbev/msn067"),)),
    ModelRow(("WAG",), "general", (Ref(
        "Whelan S, Goldman N.", "Mol Biol Evol 18(5): 691–699 (2001).",
        "10.1093/oxfordjournals.molbev.a003851"),)),
    ModelRow(("JTT",), "general", (Ref(
        "Jones DT, Taylor WR, Thornton JM.", "CABIOS 8(3): 275–282 (1992).",
        "10.1093/bioinformatics/8.3.275"),)),
    ModelRow(("JTTDCMut", "DCMut"), "general", (Ref(
        "Kosiol C, Goldman N.", "Mol Biol Evol 22(2): 193–199 (2005).",
        "10.1093/molbev/msi005"),)),
    ModelRow(("Dayhoff",), "general", (Ref(
        "Dayhoff MO, Schwartz RM, Orcutt BC.",
        "Atlas of Protein Sequence and Structure 5 (Suppl. 3): 345–352 (1978)."),)),
    ModelRow(("VT",), "general", (Ref(
        "Müller T, Vingron M.", "J Comput Biol 7(6): 761–776 (2000).",
        "10.1089/10665270050514918"),)),
    ModelRow(("PMB",), "general", (Ref(
        "Veerassamy S, Smith A, Tillier ERM.", "J Comput Biol 10(6): 997–1010 (2003).",
        "10.1089/106652703322756195"),)),
    ModelRow(("Blosum62",), "general", (Ref(
        "Henikoff S, Henikoff JG.", "PNAS 89(22): 10915–10919 (1992).",
        "10.1073/pnas.89.22.10915"),)),
    ModelRow(("Q.pfam", "Q.insect", "Q.yeast", "Q.plant", "Q.mammal", "Q.bird"),
             "estimated from large databases (Pfam; by clade)", (Ref(
        "Minh BQ, Dang CC, Vinh LS, Lanfear R.", "Syst Biol 70(5): 1046–1060 (2021).",
        "10.1093/sysbio/syab010"),)),
    ModelRow(("mtREV",), "vertebrate mitochondrial", (Ref(
        "Adachi J, Hasegawa M.", "J Mol Evol 42(4): 459–468 (1996).",
        "10.1007/BF02498640"),)),
    ModelRow(("mtART",), "arthropod mitochondrial", (Ref(
        "Abascal F, Posada D, Zardoya R.", "Mol Biol Evol 24(1): 1–5 (2007).",
        "10.1093/molbev/msl136"),)),
    ModelRow(("mtZOA",), "animal mitochondrial", (Ref(
        "Rota-Stabelli O, Yang Z, Telford MJ.",
        "Mol Phylogenet Evol 52(1): 268–272 (2009).", "10.1016/j.ympev.2009.01.011"),)),
    ModelRow(("mtMAM",), "mammalian mitochondrial", (Ref(
        "Yang Z, Nielsen R, Hasegawa M.", "Mol Biol Evol 15(12): 1600–1611 (1998).",
        "10.1093/oxfordjournals.molbev.a025888"),)),
    ModelRow(("mtInv", "mtMet", "mtVer"), "invertebrate, metazoan, vertebrate mitochondrial",
             (Ref("Le VS, Dang CC, Le QS.", "BMC Evol Biol 17: 136 (2017).",
                  "10.1186/s12862-017-0987-y"),)),
    ModelRow(("cpREV",), "chloroplast", (Ref(
        "Adachi J, Waddell PJ, Martin W, Hasegawa M.", "J Mol Evol 50(4): 348–358 (2000).",
        "10.1007/s002399910038"),)),
    ModelRow(("rtREV",), "retroviral", (Ref(
        "Dimmic MW, Rest JS, Mindell DP, Goldstein RA.", "J Mol Evol 55(1): 65–73 (2002).",
        "10.1007/s00239-001-2304-y"),)),
    ModelRow(("HIVb", "HIVw"), "HIV (between and within patients)", (Ref(
        "Nickle DC, Heath L, Jensen MA, et al.", "PLoS ONE 2(6): e503 (2007).",
        "10.1371/journal.pone.0000503"),)),
    ModelRow(("FLU",), "influenza", (Ref(
        "Dang CC, Le QS, Gascuel O, Le VS.", "BMC Evol Biol 10: 99 (2010).",
        "10.1186/1471-2148-10-99"),)),
    ModelRow(("FLAVI",), "flavivirus", (Ref(
        "Le TK, Vinh LS.", "J Mol Evol 88(5): 445–452 (2020).",
        "10.1007/s00239-020-09943-3"),)),
)  # fmt: skip

VARIANTS: tuple[ModelRow, ...] = (
    ModelRow(("+I",), "a proportion of invariant sites", (Ref(
        "Reeves JH.", "J Mol Evol 35(1): 17–31 (1992).", "10.1007/BF00160257"),)),
    ModelRow(("+G4",), "gamma-distributed rates, 4 categories", (Ref(
        "Yang Z.", "J Mol Evol 39(3): 306–314 (1994).", "10.1007/BF00160154"),)),
    ModelRow(("+I+G4",), "invariant sites plus gamma", (Ref(
        "Gu X, Fu YX, Li WH.", "Mol Biol Evol 12(4): 546–557 (1995).",
        "10.1093/oxfordjournals.molbev.a040235"),)),
    ModelRow(("+R4", "+I+R4"), "FreeRate: free rate categories (ModelFinder tries 2–10)", (
        Ref("Yang Z.", "Genetics 139(2): 993–1005 (1995).", "10.1093/genetics/139.2.993"),
        Ref("Soubrier J, et al.", "Mol Biol Evol 29(11): 3345–3358 (2012).",
            "10.1093/molbev/mss140"),
    )),
    ModelRow(("+F",), "amino-acid frequencies counted from the alignment, instead of the "
             "matrix's own"),
    ModelRow(("+C20+G4", "+C20+F+G4"), "site-profile mixture, weights fitted to the alignment; "
             "fixed models only, IQ-TREE only. +F adds a class of the alignment's own "
             "frequencies", (Ref(
        "Si Quang L, Gascuel O, Lartillot N.", "Bioinformatics 24(20): 2317–2323 (2008).",
        "10.1093/bioinformatics/btn445"),)),
)  # fmt: skip

RAXML_MISSING = sorted(catalog.RAXML_MISSING)


def cite(ref: Ref, label: str) -> str:
    """One citation as a link, for help text that must not drift from this page."""
    return f'<a href="{ref.link}">{label}</a>'


def strategy_help() -> str:
    """The "Filling the wells" popup on the library's Design section.

    It lives here, beside the citations, so the page and the popup can never cite
    different things; `_help.html` renders it with the site's own "?" affordance. The
    three strategies are PhyloBot's own, so the popup says what inspired each rather
    than implying a method paper, and reports what the simulated test measured.
    """
    return (
        "<span>All three pick from the <i>same</i> set of candidate sequences; they "
        "differ only in which of them get your wells.</span>"
        "<span><b>Most supported</b> takes the sequences that the most distinct choices "
        "of aligner, model and gap placement agree on, ranked by their summed BIC "
        "weight. Inspired by model averaging, where a result is weighted by how well "
        "each model fits rather than taken from one "
        f"({cite(_POSADA, 'Posada &amp; Buckley 2004')})."
        "</span>"
        "<span><b>Spanning</b> repeatedly takes the sequence furthest from everything "
        "already chosen (greedy max–min, or farthest-first: "
        f"{cite(_GONZALEZ, 'Gonzalez 1985')}), so the wells cover the disagreement "
        "between reconstructions instead of holding near-identical variants.</span>"
        "<span><b>Per ancestor</b> gives each ancestor an equal share of the wells, "
        "round robin, for a library over a lineage. An allocation rule, with no "
        "statistical claim.</span>"
        "<span>The three strategies are PhyloBot's own, not published methods; the "
        "citations above are the statistical inspiration. On PhyloBot's own simulated "
        "test (ten replicates with known ancestors) Spanning did not improve accuracy "
        "over Most supported, which is why Most supported is the default; what the "
        "wells are <i>for</i> — hedging against being wrong — is a different question "
        "from accuracy, so the others stay available. The sequences they choose among "
        f"come from the ML ancestor, AltAll ({cite(_EICK, 'Eick et al. 2017')}), "
        "posterior samples and single-site alternates "
        f"({cite(_VANANTWERP, 'VanAntwerp et al. 2022')})."
        "</span>"
    )


# The feature comparison under "What is distinctive": edit feature_comparison.csv (one row
# per feature, a column per tool; each cell a mark, then an optional note). Notes appear on
# hover only, unless written in [brackets], which also shows them beside the mark.
COMPARISON_FILE = Path(__file__).with_name("feature_comparison.csv")
MARKS = {"✓": "yes", "◐": "partial", "—": "no", "◌": "planned", "?": "unsure"}
COMPARISON_LEGEND = [
    ("✓", ("yes", "yes")),
    ("◐", ("partial", "partial or limited")),
    ("—", ("no", "no")),
    ("◌", ("planned", "in progress or planned")),
    ("?", ("unsure", "not yet verified")),
]


@dataclass(frozen=True)
class Mark:
    symbol: str
    kind: str  # a value of MARKS
    note: str = ""
    unsure: bool = False  # "✓?": the mark itself is not yet verified
    shown: bool = False  # the note is drawn beside the mark, not only on hover

    @property
    def title(self) -> str:
        words = [self.kind] + (["not yet verified"] if self.unsure else [])
        return ", ".join(words) + (f": {self.note}" if self.note else "")


@dataclass(frozen=True)
class Feature:
    name: str
    marks: tuple[Mark, ...]


def parse_mark(cell: str) -> Mark:
    cell = cell.strip()
    symbol = cell[:1]
    if symbol not in MARKS:
        raise ValueError(f"comparison cell {cell!r} does not start with a mark")
    rest = cell[1:]
    unsure = symbol != "?" and rest.startswith("?")
    if unsure:
        rest = rest[1:]
    note = rest.strip()
    shown = note.startswith("[") and note.endswith("]")
    if shown:
        note = note[1:-1].strip()
    return Mark(symbol, MARKS[symbol], note, unsure, shown)


type Sections = list[tuple[str, list[Feature]]]


def load_comparison(path: Path = COMPARISON_FILE) -> tuple[list[str], Sections]:
    """The tools (columns) and the features grouped by section, in the file's order."""
    with path.open(encoding="utf-8-sig", newline="") as f:
        rows = [row for row in csv.reader(f) if any(cell.strip() for cell in row)]
    tools = rows[0][2:]
    sections: Sections = []
    for section, name, *cells in rows[1:]:
        if section == "Legend":
            continue
        if len(cells) != len(tools):
            raise ValueError(f"comparison row {name!r} has {len(cells)} cells, not {len(tools)}")
        if not sections or sections[-1][0] != section:
            sections.append((section, []))
        sections[-1][1].append(Feature(name, tuple(parse_mark(c) for c in cells)))
    return tools, sections


COMPARISON_TOOLS, COMPARISON = load_comparison()
