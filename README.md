# PhyloBot 2

PhyloBot 2 automates the path from a set of protein sequences to reconstructed
ancestral proteins — alignment, phylogeny, ancestral sequence reconstruction (ASR) —
and then helps you explore those ancestors: site by site, on the tree, and along
the mutational paths that connect them.

It is a ground-up rewrite of PhyloBot (Hanson-Smith & Johnson 2016), built on
contemporary tools: Nextflow and pinned containers for the pipeline (MAFFT, MUSCLE 5,
Clustal Omega, PRANK, IQ-TREE 3, RAxML-NG, GRASP), and a small Django + HTMX web app.
It runs on a laptop, on your own servers, or in the cloud, and the pipeline also works
on its own, without the web app.

PhyloBot 2 is under active development: expect changes between releases.

## Reference

PhyloBot 2 continues the work described in:

> Hanson-Smith V, Johnson A. *PhyloBot: A Web Portal for Automated Phylogenetics,
> Ancestral Sequence Reconstruction, and Exploration of Mutational Trajectories.*
> PLoS Computational Biology 12(7): e1004976 (2016).
> [doi:10.1371/journal.pcbi.1004976](https://doi.org/10.1371/journal.pcbi.1004976)

If you use PhyloBot 2 in published work, please cite this paper, together with the
methods your analysis used (see [Included methods](#included-methods)).
The same information is in [`CITATION.cff`](CITATION.cff), which GitHub offers as
"Cite this repository" (APA and BibTeX).

## What is distinctive about PhyloBot 2

As far as we can tell from the literature and the tools we know of (see
[Related work](#related-work)), the combination below is not available elsewhere.
Several individual pieces are done well by other tools, and we say so there.

1. **Robustness across alignments, trees and models, by design.** An ancestor is an
   estimate that depends on the alignment, the tree and the substitution model behind
   it. PhyloBot 2 runs every combination you choose — up to four aligners (MAFFT,
   MUSCLE 5, Clustal Omega, PRANK) and your own alignments; tree searches by IQ-TREE
   and RAxML-NG and your own trees; any set of fixed models and a ModelFinder search —
   as separate, fully traceable analyses. It ranks the models within each alignment
   by BIC, with AIC/AICc/BIC weights, and presents the ancestors of each combination
   side by side, so you can see which conclusions hold across reasonable choices and
   which do not.

2. **A choice of gap-placement methods, each validated.** Maximum-likelihood ASR
   places a residue at every alignment column, so raw ancestors are much longer than
   any real protein. PhyloBot 2 keeps the unmasked reconstruction and adds a version
   for each method you choose: parsimony by shared presence, Fitch parsimony,
   maximum likelihood on the gap/residue pattern (with per-site gap probabilities),
   and the indel-aware reconstructions of PRANK and GRASP. Each method is scored
   against simulated protein families with known ancestors that include insertions
   and deletions, averaged over replicates, and the defaults were chosen from those
   results.

3. **Correctness checked against known answers.** The pipeline is tested
   continuously on simulated datasets whose true ancestors are known. The tests also
   rerun each dataset with the true alignment and/or true tree supplied, to separate
   error from alignment, from tree inference, and from the reconstruction itself.
   Thresholds on these tests guard every change to the pipeline.

4. **Mutational paths between ancestors.** Pick any ancestor and any of its
   descendants (another ancestor or an extant sequence), and PhyloBot 2 traces the
   path between them. It lists the substitutions, gains, and losses on each branch;
   shows every sequence along the path with each residue's posterior probability,
   outlining what changed on each branch; lets you add extant sequences for
   comparison; and exports it all to Excel. Exploring these trajectories was the
   feature the original PhyloBot introduced, and we have not found it elsewhere for
   protein ancestors. (FireProtASR 2.0 predicts *future* sequences with its Successor
   Sequence Predictor, which is a different question: what comes next, rather than what
   happened between two reconstructed nodes.)

5. **Ancestors you can inspect.**
   - **Ancestor viewer:** one ancestral node site by site, with its full posterior
     distribution at each site, summary plots of how well supported the node is
     overall and where along the sequence it is uncertain, and one or more extant
     sequences alongside.
   - **Tree viewer:** each tree with clickable ancestral nodes and branch support.
     Branch lengths can be adjusted so that very short and very long branches stay
     readable.

   Colour palettes, line styles and page skins (including dark and high-legibility
   ones) can be changed for accessibility, and plots and tables download as PNG, SVG,
   and Excel.

6. **Homologs you choose yourself.** Starting from a single protein, a BLAST search at
   EBI (UniRef90) or NCBI returns hits that you filter and review — with histograms
   of coverage and identity, and hand-picked additions or removals — before they
   become a job. Which sequences go in stays your decision.

7. **More than one kind of evidence.** For coding DNA, codon models (IQ-TREE and
   PAML codeml) run beside amino-acid models, with synonymous and non-synonymous changes
   on the mutational paths, dN/dS per branch, and a branch-site test from the tree.
   Each tree can be reconciled with a species tree (GeneRax, Open Tree of Life) to show
   duplications and losses, and GeneRax's reconciled tree can be one more tree source.
   An **ancestral library** turns ancestors and their alternatives into a short list to
   order, with AltAll controls and a methods paragraph. **Property posteriors** measure
   sequence properties (charge, hydrophobic moment, disorder, ESM-2 likelihood) over
   draws from each ancestor's posterior, so a derived number comes with an interval.
   **Simulations** evolve sequences along your own job's tree and keep the truth, so a
   job on them can be checked.

8. **Runs where you are.** The same pipeline runs from the command line with only
   Docker (and make), as a local web app (no accounts; nothing leaves your machine
   unless you choose a homolog search, which sends your query to EBI or NCBI and asks
   first), or as a hosted service. Every tool is pinned to an exact version for
   reproducibility.

## Getting started

PhyloBot 2 needs [Docker](https://docs.docker.com/get-docker/) and `make`. Every tool,
including Nextflow, runs in a container.

**The web app on your own computer.** Download PhyloBot and start it, then open
<http://localhost:8000>:

```
curl -L https://github.com/headwatersdiscovery/phylobot2/releases/latest/download/phylobot.tar.gz | tar xz
cd phylobot
make up
```

The first `make up` builds the images (several GB; allow a while), creates `.env` from
[`.env.example`](.env.example), generates this installation's own passwords and keys
in it, and sets up local storage. Its other defaults work as they are. Set
`PHYLOBOT_SEARCH_EMAIL` there if you want homolog searches (EBI and NCBI require a
contact email). Keep `.env`: your jobs' database and storage are locked with the
credentials in it.

**Upgrading.** Run `make down`, rename the old folder (`mv phylobot phylobot-old`),
download the new release with the same two lines, copy `phylobot-old/.env` into the new
`phylobot` folder, and run `make up` there. Your jobs live in Docker
volumes, not in the folder, so they carry over; the old `.env` is what lets the new
version open them.

**Check the CPU setting on every new machine.** `make docker-cpus` reads the CPUs Docker
runs on and prints the value of `PHYLOBOT_DOCKER_CPUSET` (in `.env`) that suits them;
`make up` and `make pipeline` warn when it looks wrong. Run it again after changing your
Docker runtime's CPU count. Details are in [`.env.example`](.env.example).

**On a Mac.** Install OrbStack or Docker Desktop. The bioinformatics tools are built for
amd64 only, so on Apple Silicon they run emulated, and more slowly than on an x86
machine. OrbStack's Linux VM sets aside its last CPU, and emulated IQ-TREE gets stuck
there, so every IQ-TREE step shares one core unless `PHYLOBOT_DOCKER_CPUSET` leaves it
out (`0-6` with OrbStack's 8 CPUs; `make docker-cpus` prints the value). On Docker
Desktop it is usually empty.

**On Linux** (for example an Ubuntu server):

- Install Docker Engine from [Docker's own repository](https://docs.docker.com/engine/install/ubuntu/),
  with its `docker-buildx-plugin` and `docker-compose-plugin`, not Ubuntu's `docker.io`
  package or the snap: the images need BuildKit and the Makefile calls `docker compose`.
  Use rootful Docker: the stack mounts `/var/run/docker.sock`.
- Add yourself to the `docker` group so `make` runs without `sudo`:
  `sudo usermod -aG docker $USER`, then log out and back in (reconnect over SSH; start
  new tmux or screen sessions). `id -nG` should list `docker`, and `docker ps` should
  work without `sudo`.
- Leave `PHYLOBOT_DOCKER_CPUSET` empty: Docker runs on your own kernel, which sets no CPU
  aside. Set it only to keep pipeline steps off some cores, such as `0-13` on a
  16-thread machine you also use for other work. On x86 the tools run natively, with
  no emulation.
- Pipeline containers run as root and write into the project folder (`jobs/`, `work/`,
  `results/`, `.nextflow*`), so those files belong to root: use `sudo make clean`.
- Postgres and Garage keep their data in Docker volumes (under `/var/lib/docker`), and
  the first `make up` builds several GB of images: put that on a disk with room.
- **Using it from another computer.** Local mode has no accounts and listens only on
  the machine itself (127.0.0.1); it refuses other hostnames. Do not expose it on your
  network. From your laptop, open an SSH tunnel and browse as if it were local:

  ```
  ssh -L 8000:127.0.0.1:8000 you@server
  ```

  then open <http://localhost:8000> on the laptop.

Local mode has no accounts and listens only on your own machine. **New job** lets you
paste or upload sequences (or load one of the bundled example families), choose aligners, tree
searches, models, an outgroup and gap-placement methods, and see an estimate of the
work before you submit.

**The pipeline alone, without the web app.** From an unaligned protein FASTA file:

```
make pipeline IN=my_proteins.fasta OUT=results/my_proteins \
  ARGS="--aligners mafft,muscle --tree_builders iqtree,raxml --outgroup Some.outgroup.seq"
```

From coding DNA, add `--seqtype codon` (and `--genetic_code`, an NCBI number; 1, the
standard code, by default): the sequences are checked (whole codons, no internal stop
codons) and translated, the protein steps run on the translations, and each aligner's
codon alignment is published as `<aligner>/alignment.cds.fasta`. The `prank_codon`
aligner aligns the DNA directly with PRANK's codon model. `--codon_models` (for
example `GY+F3X4+G4`; the list is in `src/phylobot/catalog.py`) adds combinations that
run IQ-TREE under codon models on the codon alignment: a tree search, fits to your own
trees, and ancestors as codons (`ancestors.cds.fasta`) and as proteins.

**Simulations.** The Simulations page (top bar) evolves sequences along one of your
jobs' trees or a random tree and keeps the truth; each replicate can be used in a new
job, whose results can then be checked against it. From the command line,
`--simulate spec.json` simulates sequence evolution with AliSim and
writes each replicate's truth (unaligned input, true alignment, ancestors and tree,
provenance) to `<OUT>/replicates/<seed>/`; a spec on your own tree or root sequence takes
`--sim_tree` / `--sim_root`. Spec format: `src/phylobot/simulate/spec.py` (the golden
datasets' `tests/golden/*/*/spec.toml` are examples, in TOML). A spec with a `[family]`
section reads `--sim_tree` as a rooted species tree and first grows a gene family on it
with duplications and losses (random rates and/or chosen events), then evolves sequences
along the gene tree; the true events and each gene's species join the replicate's truth.

Results go to `<OUT>/<aligner>/<tree>/<model>/` (`ancestors.fasta`, per-site
posterior probabilities, `tree.nwk`, model fit statistics, and one
`ancestors.gaps-<method>.fasta` per gap-placement method), with `summary.json` and
`model_comparison.tsv` at the top. Every option is documented in
[`pipeline/nextflow.config`](pipeline/nextflow.config).

**Developers:** `git clone https://github.com/headwatersdiscovery/phylobot2.git`. You are welcome
to edit or extend PhyloBot's code under its [licence](#licence). See
[Development](#development).

## Included methods

PhyloBot 2 is a workflow around other people's methods. Please cite the ones your
analysis used. Versions are those pinned in
[`pipeline/nextflow.config`](pipeline/nextflow.config).

### Multiple sequence alignment

- **MAFFT** 7.525 — Katoh K, Standley DM. *MAFFT Multiple Sequence Alignment Software
  Version 7: Improvements in Performance and Usability.* Molecular Biology and
  Evolution 30(4): 772–780 (2013).
  [doi:10.1093/molbev/mst010](https://doi.org/10.1093/molbev/mst010)
- **MUSCLE** 5.3 — Edgar RC. *Muscle5: High-accuracy alignment ensembles enable
  unbiased assessments of sequence homology and phylogeny.* Nature Communications 13:
  6968 (2022). [doi:10.1038/s41467-022-34630-w](https://doi.org/10.1038/s41467-022-34630-w)
- **Clustal Omega** 1.2.4 — Sievers F, Wilm A, Dineen D, et al. *Fast, scalable
  generation of high-quality protein multiple sequence alignments using Clustal
  Omega.* Molecular Systems Biology 7: 539 (2011).
  [doi:10.1038/msb.2011.75](https://doi.org/10.1038/msb.2011.75)
- **PRANK** (v.250331) — Löytynoja A, Goldman N. *An algorithm for progressive
  multiple alignment of sequences with insertions.* PNAS 102(30): 10557–10562 (2005).
  [doi:10.1073/pnas.0409137102](https://doi.org/10.1073/pnas.0409137102); and
  Löytynoja A, Goldman N. *Phylogeny-aware gap placement prevents errors in sequence
  alignment and evolutionary analysis.* Science 320(5883): 1632–1635 (2008).
  [doi:10.1126/science.1158395](https://doi.org/10.1126/science.1158395)

### Phylogeny, model selection and branch support

- **IQ-TREE** 3.1.3 (tree search, model fitting, ancestral reconstruction) — Wong TKF,
  Ly-Trong N, Ren H, et al. *IQ-TREE 3: phylogenomic inference software using complex
  evolutionary models.* Molecular Biology and Evolution 43(5): msag117 (2026).
  [doi:10.1093/molbev/msag117](https://doi.org/10.1093/molbev/msag117)
- **RAxML-NG** 2.0.3 (a second tree search; its trees are then fitted by IQ-TREE) —
  Kozlov AM, Darriba D, Flouri T, Morel B, Stamatakis A. *RAxML-NG: a fast, scalable
  and user-friendly tool for maximum likelihood phylogenetic inference.*
  Bioinformatics 35(21): 4453–4455 (2019).
  [doi:10.1093/bioinformatics/btz305](https://doi.org/10.1093/bioinformatics/btz305);
  and Kozlov O, Togkousidis A, Stelz C, et al. *RAxML-NG 2: Automatic model
  selection, novel tree search heuristics, and fast branch support metrics.* bioRxiv
  preprint (2026). [doi:10.64898/2026.09.09.750097](https://doi.org/10.64898/2026.09.09.750097)
- **GeneRax** 2.1.3 (reconciles each tree with the species tree; searches a reconciled
  tree) — Morel B, Kozlov AM, Stamatakis A, Szöllősi GJ. *GeneRax: a tool for
  species-tree-aware maximum likelihood-based gene family tree inference under gene
  duplication, transfer, and loss.* Molecular Biology and Evolution 37(9): 2763–2774
  (2020). [doi:10.1093/molbev/msaa141](https://doi.org/10.1093/molbev/msaa141)
- **Open Tree of Life** (the species tree for reconciliation: synthetic tree and
  taxonomy; the release is recorded with each job) — Hinchliff CE, Smith SA, Allman JF,
  et al. *Synthesis of phylogeny and taxonomy into a comprehensive tree of life.*
  Proceedings of the National Academy of Sciences 112(41): 12764–12769 (2015).
  [doi:10.1073/pnas.1423041112](https://doi.org/10.1073/pnas.1423041112); and
  Redelings BD, Holder MT. *A supertree pipeline for summarizing phylogenetic and
  taxonomic information for millions of species.* PeerJ 5: e3058 (2017).
  [doi:10.7717/peerj.3058](https://doi.org/10.7717/peerj.3058)
- **UniProt** (species read from UniProt entry names, using its controlled vocabulary of
  species, CC BY 4.0, bundled; and from UniProt accessions, with consent) — The UniProt
  Consortium. *UniProt: the Universal Protein Knowledgebase in 2025.* Nucleic Acids
  Research 53(D1): D609–D617 (2025).
  [doi:10.1093/nar/gkae1010](https://doi.org/10.1093/nar/gkae1010)
- **AU test** (in IQ-TREE; is each tree source's topology rejected by the sequences?) —
  Shimodaira H. *An approximately unbiased test of phylogenetic tree selection.*
  Systematic Biology 51(3): 492–508 (2002).
  [doi:10.1080/10635150290069913](https://doi.org/10.1080/10635150290069913)
- **ModelFinder** (in IQ-TREE) — Kalyaanamoorthy S, Minh BQ, Wong TKF, von Haeseler A,
  Jermiin LS. *ModelFinder: fast model selection for accurate phylogenetic
  estimates.* Nature Methods 14: 587–589 (2017).
  [doi:10.1038/nmeth.4285](https://doi.org/10.1038/nmeth.4285)
- **Model comparison** by AIC — Akaike H, IEEE Transactions on Automatic Control
  19(6): 716–723 (1974), [doi:10.1109/TAC.1974.1100705](https://doi.org/10.1109/TAC.1974.1100705);
  AICc — Hurvich CM, Tsai CL, Biometrika 76(2): 297–307 (1989),
  [doi:10.1093/biomet/76.2.297](https://doi.org/10.1093/biomet/76.2.297);
  and BIC — Schwarz G, Annals of Statistics 6(2): 461–464 (1978),
  [doi:10.1214/aos/1176344136](https://doi.org/10.1214/aos/1176344136).
- **UFBoot** (ultrafast bootstrap) — Hoang DT, Chernomor O, von Haeseler A, Minh BQ,
  Vinh LS. *UFBoot2: Improving the Ultrafast Bootstrap Approximation.* Molecular
  Biology and Evolution 35(2): 518–522 (2018).
  [doi:10.1093/molbev/msx281](https://doi.org/10.1093/molbev/msx281)
- **SH-aLRT** (SH-like approximate likelihood ratio test) — Guindon S, Dufayard JF,
  Lefort V, Anisimova M, Hordijk W, Gascuel O. *New Algorithms and Methods to Estimate
  Maximum-Likelihood Phylogenies: Assessing the Performance of PhyML 3.0.* Systematic
  Biology 59(3): 307–321 (2010).
  [doi:10.1093/sysbio/syq010](https://doi.org/10.1093/sysbio/syq010)

### Ancestral sequence reconstruction

- **Empirical Bayes (marginal) reconstruction**, as implemented in IQ-TREE
  (`--ancestral`) — Yang Z, Kumar S, Nei M. *A new method of inference of ancestral
  nucleotide and amino acid sequences.* Genetics 141(4): 1641–1650 (1995).
  [doi:10.1093/genetics/141.4.1641](https://doi.org/10.1093/genetics/141.4.1641)

### Gaps in ancestral sequences

- **Parsimony (shared presence)** — PhyloBot 2's own rule: a residue only where it is
  present on both sides of a node.
- **Parsimony (Fitch)** — Fitch WM. *Toward Defining the Course of Evolution: Minimum
  Change for a Specific Tree Topology.* Systematic Zoology 20(4): 406–416 (1971).
  [doi:10.2307/2412116](https://doi.org/10.2307/2412116); multifurcations by
  Hartigan JA. *Minimum Mutation Fits to a Given Tree.* Biometrics 29(1): 53–65
  (1973). [doi:10.2307/2529676](https://doi.org/10.2307/2529676)
- **Maximum likelihood** — the gap/residue pattern reconstructed by IQ-TREE under a
  two-state model on the fixed tree (IQ-TREE 3, above).
- **PRANK** ancestors (`-showanc`) — PRANK, above.
- **GRASP** 2.0.0 (joint reconstruction with bi-directional edge parsimony for indels)
  — Foley G, Mora A, Ross CM, et al. *Engineering indel and substitution variants of
  diverse and ancient enzymes using Graphical Representation of Ancestral Sequence
  Predictions (GRASP).* PLoS Computational Biology 18(10): e1010633 (2022).
  [doi:10.1371/journal.pcbi.1010633](https://doi.org/10.1371/journal.pcbi.1010633)

### Sequence properties

Measured over draws from each ancestor's posterior.

- **Physicochemical descriptors** (net charge, hydrophobic moment, hydrophobic regions),
  summarised by PhyloBot's own windows and thresholds. Hydrophobic moment — Eisenberg D, Weiss RM, Terwilliger TC. *The helical
  hydrophobic moment: a measure of the amphiphilicity of a helix.* Nature 299: 371–374
  (1982). [doi:10.1038/299371a0](https://doi.org/10.1038/299371a0); on the consensus
  scale of Eisenberg D, Schwarz E, Komaromy M, Wall R. Journal of Molecular Biology
  179(1): 125–142 (1984).
  [doi:10.1016/0022-2836(84)90309-7](https://doi.org/10.1016/0022-2836(84)90309-7)
- **Hydrophobic regions** — Kyte J, Doolittle RF. *A simple method for displaying the
  hydropathic character of a protein.* Journal of Molecular Biology 157(1): 105–132
  (1982). [doi:10.1016/0022-2836(82)90515-0](https://doi.org/10.1016/0022-2836(82)90515-0)
- **metapredict** 3.0.2 (disordered fraction, optional) — Emenecker RJ, Griffith D,
  Holehouse AS. *Metapredict: a fast, accurate, and easy-to-use predictor of consensus
  disorder and structure.* Biophysical Journal 120(20): 4312–4319 (2021).
  [doi:10.1016/j.bpj.2021.08.039](https://doi.org/10.1016/j.bpj.2021.08.039); its V3
  network: Lotthammer JM, Hernández-García J, Griffith D, Weijers D, Holehouse AS,
  Emenecker RJ. *Metapredict enables accurate disorder prediction across the Tree of
  Life.* bioRxiv (2024), preprint.
  [doi:10.1101/2024.11.05.622168](https://doi.org/10.1101/2024.11.05.622168)
- **ESM-2** 35M (log-likelihood, optional) — Lin Z, Akin H, Rao R, et al.
  *Evolutionary-scale prediction of atomic-level protein structure with a language
  model.* Science 379(6637): 1123–1130 (2023).
  [doi:10.1126/science.ade2574](https://doi.org/10.1126/science.ade2574)

### Homolog search

- **BLAST** (blastp) — Camacho C, Coulouris G, Avagyan V, et al. *BLAST+: architecture
  and applications.* BMC Bioinformatics 10: 421 (2009).
  [doi:10.1186/1471-2105-10-421](https://doi.org/10.1186/1471-2105-10-421)
- **EMBL-EBI Job Dispatcher** — Madeira F, Madhusoodanan N, Lee J, et al. *The
  EMBL-EBI Job Dispatcher sequence analysis tools framework in 2024.* Nucleic Acids
  Research 52(W1): W521–W525 (2024).
  [doi:10.1093/nar/gkae241](https://doi.org/10.1093/nar/gkae241)
- **UniRef90** — Suzek BE, Wang Y, Huang H, McGarvey PB, Wu CH, UniProt Consortium.
  *UniRef clusters: a comprehensive and scalable alternative for improving sequence
  similarity searches.* Bioinformatics 31(6): 926–932 (2015).
  [doi:10.1093/bioinformatics/btu739](https://doi.org/10.1093/bioinformatics/btu739)
- **NCBI BLAST** with ClusteredNR (local mode only) —
  [blast.ncbi.nlm.nih.gov](https://blast.ncbi.nlm.nih.gov/)

### Workflow and testing

- **Nextflow** — Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C.
  *Nextflow enables reproducible computational workflows.* Nature Biotechnology 35:
  316–319 (2017). [doi:10.1038/nbt.3820](https://doi.org/10.1038/nbt.3820)
- **BioContainers** (the pinned tool images) — da Veiga Leprevost F, Grüning BA, Alves
  Aflitos S, et al. *BioContainers: an open-source and community-driven framework for
  software standardization.* Bioinformatics 33(16): 2580–2582 (2017).
  [doi:10.1093/bioinformatics/btx192](https://doi.org/10.1093/bioinformatics/btx192)
- **DendroPy** (tree handling in PhyloBot's own steps) — Moreno MA, Holder MT,
  Sukumaran J. *DendroPy 5: a mature Python library for phylogenetic computing.*
  Journal of Open Source Software 9(101): 6943 (2024).
  [doi:10.21105/joss.06943](https://doi.org/10.21105/joss.06943)
- **AliSim** (the simulated test datasets with known ancestors) — Ly-Trong N,
  Naser-Khdour S, Lanfear R, Minh BQ. *AliSim: A Fast and Versatile Phylogenetic
  Sequence Simulator for the Genomic Era.* Molecular Biology and Evolution 39(5):
  msac092 (2022). [doi:10.1093/molbev/msac092](https://doi.org/10.1093/molbev/msac092)

### Substitution models

Models are chosen from a fixed catalog, as an amino-acid matrix plus a rate and
frequency variant (for example `LG+F+G4`). Every one can be used as a fixed model or
as a ModelFinder candidate; RAxML-NG lacks JTTDCMut.

| Matrix | Use | Reference |
|---|---|---|
| LG | general | Le SQ, Gascuel O. MBE 25(7): 1307–1320 (2008). [doi:10.1093/molbev/msn067](https://doi.org/10.1093/molbev/msn067) |
| WAG | general | Whelan S, Goldman N. MBE 18(5): 691–699 (2001). [doi:10.1093/oxfordjournals.molbev.a003851](https://doi.org/10.1093/oxfordjournals.molbev.a003851) |
| JTT | general | Jones DT, Taylor WR, Thornton JM. CABIOS 8(3): 275–282 (1992). [doi:10.1093/bioinformatics/8.3.275](https://doi.org/10.1093/bioinformatics/8.3.275) |
| JTTDCMut, DCMut | general | Kosiol C, Goldman N. MBE 22(2): 193–199 (2005). [doi:10.1093/molbev/msi005](https://doi.org/10.1093/molbev/msi005) |
| Dayhoff | general | Dayhoff MO, Schwartz RM, Orcutt BC. In: Atlas of Protein Sequence and Structure 5 (Suppl. 3): 345–352 (1978) |
| VT | general | Müller T, Vingron M. J Comput Biol 7(6): 761–776 (2000). [doi:10.1089/10665270050514918](https://doi.org/10.1089/10665270050514918) |
| PMB | general | Veerassamy S, Smith A, Tillier ERM. J Comput Biol 10(6): 997–1010 (2003). [doi:10.1089/106652703322756195](https://doi.org/10.1089/106652703322756195) |
| Blosum62 | general | Henikoff S, Henikoff JG. PNAS 89(22): 10915–10919 (1992). [doi:10.1073/pnas.89.22.10915](https://doi.org/10.1073/pnas.89.22.10915) |
| Q.pfam, Q.insect, Q.yeast, Q.plant, Q.mammal, Q.bird | estimated from large databases (Pfam; by clade) | Minh BQ, Dang CC, Vinh LS, Lanfear R. Syst Biol 70(5): 1046–1060 (2021). [doi:10.1093/sysbio/syab010](https://doi.org/10.1093/sysbio/syab010) |
| mtREV | vertebrate mitochondrial | Adachi J, Hasegawa M. J Mol Evol 42(4): 459–468 (1996). [doi:10.1007/BF02498640](https://doi.org/10.1007/BF02498640) |
| mtART | arthropod mitochondrial | Abascal F, Posada D, Zardoya R. MBE 24(1): 1–5 (2007). [doi:10.1093/molbev/msl136](https://doi.org/10.1093/molbev/msl136) |
| mtZOA | animal mitochondrial | Rota-Stabelli O, Yang Z, Telford MJ. Mol Phylogenet Evol 52(1): 268–272 (2009). [doi:10.1016/j.ympev.2009.01.011](https://doi.org/10.1016/j.ympev.2009.01.011) |
| mtMAM | mammalian mitochondrial | Yang Z, Nielsen R, Hasegawa M. MBE 15(12): 1600–1611 (1998). [doi:10.1093/oxfordjournals.molbev.a025888](https://doi.org/10.1093/oxfordjournals.molbev.a025888) |
| mtInv, mtMet, mtVer | invertebrate, metazoan, vertebrate mitochondrial | Le VS, Dang CC, Le QS. BMC Evol Biol 17: 136 (2017). [doi:10.1186/s12862-017-0987-y](https://doi.org/10.1186/s12862-017-0987-y) |
| cpREV | chloroplast | Adachi J, Waddell PJ, Martin W, Hasegawa M. J Mol Evol 50(4): 348–358 (2000). [doi:10.1007/s002399910038](https://doi.org/10.1007/s002399910038) |
| rtREV | retroviral | Dimmic MW, Rest JS, Mindell DP, Goldstein RA. J Mol Evol 55(1): 65–73 (2002). [doi:10.1007/s00239-001-2304-y](https://doi.org/10.1007/s00239-001-2304-y) |
| HIVb, HIVw | HIV (between and within patients) | Nickle DC, Heath L, Jensen MA, et al. PLoS ONE 2(6): e503 (2007). [doi:10.1371/journal.pone.0000503](https://doi.org/10.1371/journal.pone.0000503) |
| FLU | influenza | Dang CC, Le QS, Gascuel O, Le VS. BMC Evol Biol 10: 99 (2010). [doi:10.1186/1471-2148-10-99](https://doi.org/10.1186/1471-2148-10-99) |
| FLAVI | flavivirus | Le TK, Vinh LS. J Mol Evol 88(5): 445–452 (2020). [doi:10.1007/s00239-020-09943-3](https://doi.org/10.1007/s00239-020-09943-3) |

| Variant | Meaning | Reference |
|---|---|---|
| +I | a proportion of invariant sites | Reeves JH. J Mol Evol 35(1): 17–31 (1992). [doi:10.1007/BF00160257](https://doi.org/10.1007/BF00160257) |
| +G4 | gamma-distributed rates, 4 categories | Yang Z. J Mol Evol 39(3): 306–314 (1994). [doi:10.1007/BF00160154](https://doi.org/10.1007/BF00160154) |
| +I+G4 | invariant sites plus gamma | Gu X, Fu YX, Li WH. MBE 12(4): 546–557 (1995). [doi:10.1093/oxfordjournals.molbev.a040235](https://doi.org/10.1093/oxfordjournals.molbev.a040235) |
| +R4, +I+R4 | FreeRate: free rate categories (ModelFinder tries 2–10) | Yang Z. Genetics 139(2): 993–1005 (1995), [doi:10.1093/genetics/139.2.993](https://doi.org/10.1093/genetics/139.2.993); Soubrier J, et al. MBE 29(11): 3345–3358 (2012), [doi:10.1093/molbev/mss140](https://doi.org/10.1093/molbev/mss140) |
| +F | amino-acid frequencies counted from the alignment, instead of the matrix's own | — |

## Related work

Ancestral sequence reconstruction has an active community, and PhyloBot 2 stands on
its work. These are the tools and methods closest to what PhyloBot 2 does, with an
honest account of how they compare. Several are ahead of PhyloBot 2 in important
respects.

**PhyloBot** — Hanson-Smith V, Johnson A. *PhyloBot: A Web Portal for Automated
Phylogenetics, Ancestral Sequence Reconstruction, and Exploration of Mutational
Trajectories.* PLoS Computational Biology 12(7): e1004976 (2016).
[doi:10.1371/journal.pcbi.1004976](https://doi.org/10.1371/journal.pcbi.1004976)
The original. It automated alignment with several aligners, model fitting, branch
support, and ASR, and introduced interactive exploration of mutational trajectories.
PhyloBot 2 keeps those goals and replaces the implementation. The original's
prediction of functional loci has not yet been rebuilt.

**FireProtASR** — Musil M, Khan RT, Beier A, Stourac J, Konegger H, Damborsky J,
Bednar D. *FireProtASR: A Web Server for Fully Automated Ancestral Sequence
Reconstruction.* Briefings in Bioinformatics 22(4): bbaa337 (2021).
[doi:10.1093/bib/bbaa337](https://doi.org/10.1093/bib/bbaa337) — and **FireProtASR
2.0**, Kohout P, Lacko D, Musil M, et al., Briefings in Bioinformatics 27(4): bbag436
(2026). [doi:10.1093/bib/bbag436](https://doi.org/10.1093/bib/bbag436)
The closest in spirit: a one-stop web server that can start from a single protein
sequence, find homologs, align them (Clustal Omega), build and root the tree
(RAxML; rooting without an outgroup), reconstruct ancestors (Lazarus/PAML), and
place gaps with its own algorithm, showing per-site probabilities as logos. It is
excellent for users who want an answer from one sequence with no expertise, and its
essential-residue filter — keeping only homologs that carry a protein's catalytic or
binding positions — is a good idea we intend to generalise. Version 2.0 adds two
methods that look *forward* rather than back: the Successor Sequence Predictor, which
extrapolates physicochemical trends along lineages to predict a protein's next step,
and a variational autoencoder that samples ancestor-like sequences from an
evolutionary latent space. Both are on our roadmap.
PhyloBot 2 can start from one sequence too, but not automatically: a BLAST search at
EBI (UniRef90) or NCBI returns hits that you filter and review (with histograms of
coverage and identity, and hand-picked additions or removals) before they become a
job; which sequences go in stays your decision.
In turn, FireProtASR uses one aligner and one pipeline configuration per job, and does
not compare alignments or models, trace the path between two reconstructed ancestors,
or run outside its web server.

**EnzymeMiner** — Hon J, Borko S, Stourac J, et al. *EnzymeMiner: automated mining of
soluble enzymes with diverse structures, catalytic properties and stabilities.* Nucleic
Acids Research 48: W104–W109 (2020).
[doi:10.1093/nar/gkaa372](https://doi.org/10.1093/nar/gkaa372) — and **EnzymeMiner
2.0**, Rosinska M, Svobodova L, Borko S, et al., Nucleic Acids Research 54: W257–W265
(2026). [doi:10.1093/nar/gkag424](https://doi.org/10.1093/nar/gkag424)
Not an ASR tool, but the engine behind FireProtASR's sequence collection, and
**complementary to PhyloBot 2 rather than competing with it**. EnzymeMiner mines
enzyme sequence space exhaustively — version 2.0 searches NCBI nr (~850 million) plus
MGnify (2.4 billion) over 6 to 36 hours, annotating hits with predicted melting
temperature, solubility, optimal pH and catalytic efficiency. PhyloBot 2's homolog
search is deliberately the opposite: one BLAST search of UniRef90 at EBI (up to 1,000
hits, in minutes to a quarter of an hour), chosen for clean taxonomy and an
interpretable tree rather than breadth. Anyone needing exhaustive
mining should run EnzymeMiner and bring the resulting FASTA to PhyloBot, which already
accepts one. Their two-stage validation of essential residues — pairwise, then
re-checked within a multiple alignment — is the right way to do that filter, and we
plan to follow it.

**Topiary** — Orlandi KN, Phillips SR, Sailer ZR, Harman JL, Harms MJ. *Topiary:
Pruning the manual labor from ancestral sequence reconstruction.* Protein Science
32(2): e4551 (2023). [doi:10.1002/pro.4551](https://doi.org/10.1002/pro.4551)
A Python and command-line pipeline that excels at the hardest manual step:
building a good dataset. It finds sequences by BLAST, does taxonomically informed
quality control and redundancy reduction, reconciles the gene tree with the species
tree (GeneRax), and reconstructs ancestors (Muscle 5, RAxML-NG, PastML), producing
an HTML summary and spreadsheets. PhyloBot 2 has nothing comparable for dataset
construction or tree reconciliation. Topiary uses one aligner, has no web interface,
and does not offer a choice of gap-placement methods or trajectory exploration.

**GRASP** — Foley G, et al. *Engineering indel and substitution variants of diverse
and ancient enzymes using Graphical Representation of Ancestral Sequence Predictions
(GRASP).* PLoS Computational Biology 18(10): e1010633 (2022).
[doi:10.1371/journal.pcbi.1010633](https://doi.org/10.1371/journal.pcbi.1010633)
The leading tool for insertions and deletions. GRASP represents ancestors as partial
order graphs, infers indels explicitly, scales to families of more than 10,000
sequences, suggests indel variants for protein engineering, and has a web
interface for exploring ancestors. It starts from your alignment and tree rather
than building them. PhyloBot 2 runs GRASP (2.0.0) as one of its gap-placement
methods; in our simulations its ancestors are close to the most accurate we
tested. We do not yet
expose GRASP's alternative ancestors or its engineering features.

**AP-LASR** — VanAntwerp J, Mardikoraem M, Pascual N, Woldring D. *AP-LASR:
Automated Protein Libraries from Ancestral Sequence Reconstruction.* bioRxiv preprint
2023.10.09.561537 (2023).
[doi:10.1101/2023.10.09.561537](https://doi.org/10.1101/2023.10.09.561537)
The closest published work to PhyloBot 2's **ancestral library**: it automates ASR
specifically in order to generate libraries of ancestral and near-ancestral sequences,
with the ancestral probability cut-off as the control over how diverse the library is.
The same group's earlier chapter on alternate amino acid states for library design
(VanAntwerp et al., Methods Mol Biol 2491: 75–86, 2022,
[doi:10.1007/978-1-0716-2285-8_4](https://doi.org/10.1007/978-1-0716-2285-8_4)) is the
published precedent for PhyloBot 2's single-site alternates axis.

**AncFlow** — *AncFlow: An Ancestral Sequence Reconstruction Approach for
Determining Novel Protein Structural Landscapes.* bioRxiv preprint (2024).
[doi:10.1101/2024.07.30.605920](https://doi.org/10.1101/2024.07.30.605920)
A Snakemake pipeline (MAFFT, IQ-TREE, subfamily clustering with AutoPhy, ASR) that
reconstructs the ancestors of automatically identified subfamilies. Like PhyloBot 2,
it builds on MAFFT and IQ-TREE in a workflow manager; it is a batch pipeline rather
than an interactive tool.

**ARPIP** — Jowkar G, Pečerska J, Maiolo M, Gil M, Anisimova M. *ARPIP: Ancestral
Sequence Reconstruction with Insertions and Deletions under the Poisson Indel
Process.* Systematic Biology 72(2): 307–318 (2023).
[doi:10.1093/sysbio/syac050](https://doi.org/10.1093/sysbio/syac050)
**indelMaP** — Iglhaut C, Pečerska J, Gil M, Anisimova M. *Please Mind the Gap:
Indel-Aware Parsimony for Fast and Accurate Ancestral Sequence Reconstruction and
Multiple Sequence Alignment Including Long Indels.* Molecular Biology and Evolution
41(7): msae109 (2024). [doi:10.1093/molbev/msae109](https://doi.org/10.1093/molbev/msae109)
Reconstruction methods with explicit models of insertion and deletion, which place
indel events on the tree in a principled way. They are methods rather than end-user
tools, and strong candidates to add as gap-placement methods in PhyloBot 2, where
our simulations could compare them with the others.

**Nextstrain / Auspice** — Hadfield J, et al. *Nextstrain: real-time tracking of
pathogen evolution.* Bioinformatics 34(23): 4121–4123 (2018).
[doi:10.1093/bioinformatics/bty407](https://doi.org/10.1093/bioinformatics/bty407)
An outstanding interactive tree viewer that shows the nucleotide and amino-acid
mutations on each branch, from ancestral reconstruction. It is built for tracking
pathogens across many genomes, not for protein ancestors: it does not show
posterior probabilities or compare reconstruction settings. Its branch-mutation
display is the nearest relative of PhyloBot 2's mutational paths.

**MEGA 12** — Kumar S, Stecher G, Suleski M, Sanderford M, Sharma S, Tamura K.
*MEGA12: Molecular Evolutionary Genetic Analysis Version 12 for Adaptive and Green
Computing.* Molecular Biology and Evolution 41(12): msae263 (2024).
[doi:10.1093/molbev/msae263](https://doi.org/10.1093/molbev/msae263)
A widely used desktop application that includes maximum-likelihood ancestral states
within a broad evolutionary-analysis toolkit. It is general-purpose and interactive,
but not an automated ASR pipeline and not built for comparing reconstructions.

**Recent research directions.** Much recent work improves the reconstruction itself
rather than the workflow around it. Examples include the robustness of ASR to
evolutionary heterogeneity and epistasis — Muñiz-Trejo R, Park Y, Thornton JW,
Molecular Biology and Evolution 42(4): msaf084 (2025),
[doi:10.1093/molbev/msaf084](https://doi.org/10.1093/molbev/msaf084) — and
generative-model and coevolution-aware approaches to ASR (2026 preprints on
[bioRxiv](https://www.biorxiv.org/content/10.64898/2026.01.18.700141.full.pdf) and
[arXiv](https://arxiv.org/pdf/2606.27942)). PhyloBot 2 aims to make methods like
these easy to run, compare, and inspect as they mature.

This survey reflects the literature and tools we know of as of September 2026 and is
surely incomplete. If we have missed or misdescribed your work, please contact us at [hello@phylobot.com](mailto:hello@phylobot.com).

## Development

- Bug reports, questions and pull requests are welcome on
  [GitHub](https://github.com/headwatersdiscovery/phylobot2/issues).
- The web app (`src/phylobot/`, Django) only queues jobs and shows results; a worker runs
  each job's Nextflow pipeline (`pipeline/main.nf`), every step in a pinned container,
  and keeps sequences, trees and results in object storage, never in the database.
- `make test` runs the fast tests in a container; `make golden` runs the pipeline on
  the simulated datasets in `tests/golden/` and checks its accuracy; `make e2e` drives
  the running web stack end to end; `make lint` and `make typecheck` check the code.
- Python 3.14 with [uv](https://docs.astral.sh/uv/); Django 6 + HTMX; Postgres;
  S3-compatible object storage (Garage locally).

## Licence

PhyloBot 2 is released under the [BSD 3-Clause License](LICENSE). The tools it runs
keep their own licences; they run as separate programs in their own containers.
GRASP's code has no licence file, so PhyloBot builds its GRASP image locally and never
redistributes it.
