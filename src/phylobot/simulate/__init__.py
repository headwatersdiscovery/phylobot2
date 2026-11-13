"""Simulating sequence evolution with known truth.

One engine for every simulation PhyloBot runs: the golden datasets
(`tests/golden/generate.py`), the work-estimate benchmarks and users' own simulations.
A simulation is described by a `Spec` (`spec.py`: tree, root sequence, model, indels,
seeds), validated before anything runs, and is run in three steps so that AliSim always
runs in its own pinned container:

    prepare  (Python)  writes AliSim's inputs and a plan of its commands (plan.json)
    AliSim   (IQ-TREE) runs each planned command
    collect  (Python)  reads AliSim's outputs back into a truth bundle

The pipeline runs them as SIM_PREPARE, SIM_ALISIM and SIM_COLLECT (`--simulate`); the
golden generator runs them on the host with `engine.simulate_locally`. Either way a
replicate's truth bundle is the golden layout:

    input.fasta           unaligned leaf sequences: the only thing the pipeline sees
    true_alignment.fasta  leaves as simulated, aligned
    true_ancestors.fasta  every internal node, aligned to true_alignment.fasta
    true_tree.nwk         the simulated tree, internal nodes labelled to match
    true_site_classes.tsv each codon site's dN/dS class (site-class specs only)
    names.tsv             simulated name -> original name (neutral names only)
    species_tree.nwk,     gene families only (`[family]`, families.py): the species
    species.tsv,          tree, each gene's species, every duplication, loss and
    true_events.tsv,      speciation with the gene node it became, and counts
    family.json
    provenance.json       tool version, exact command, and checksums of the above
"""

from phylobot.simulate.spec import Spec, SpecError, load, parse

__all__ = ["Spec", "SpecError", "load", "parse"]
