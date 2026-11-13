#!/usr/bin/env nextflow
/*
 * PhyloBot: protein sequences -> alignments -> ML trees -> ancestral sequences.
 *
 * --seqtype codon takes coding DNA instead (--genetic_code, NCBI number, default 1).
 * VALIDATE translates it; every amino-acid step runs on the translations unchanged, and
 * CODON_ALIGNMENT back-translates each protein alignment (or translates PRANK -codon's,
 * aligner 'prank_codon'), publishing <aligner>/alignment.cds.fasta beside
 * alignment.fasta. --codon_models (catalog codon models: src/phylobot/catalog.py) adds
 * combinations that run on the codon alignment (meta.codon): TREE_IQTREE searches under
 * the codon model, FIT_IQTREE fits it to provided trees, ANCESTRAL_IQTREE reconstructs
 * codons, and EXPORT sums them into amino acids (ancestors.cds.fasta beside). Gap
 * placement and every export get the protein alignment (one column per codon). RAxML-NG
 * has no codon models, so codon combinations get no RAxML-NG tree. --codon_modelfinder
 * adds 'modelfinder_codon' combinations: MODELFINDER_CODON picks among --codon_candidates
 * once per codon alignment, for every tree source (as MODELFINDER does for proteins).
 * --codon_trees chooses the codon combinations' IQ-TREE trees: 'search' (a tree search
 * under each codon model; the default) and/or 'amino' (FIT_CODON_AA: each codon model
 * fitted on an amino-acid tree, no codon search). --codon_amino_trees: 'best' (the
 * default) fits on each alignment's amino-acid tree of lowest BIC per tree source
 * (tree key iqtree_aa, raxml_aa), 'all' on every amino-acid tree (iqtree_aa-<model>).
 * Provided trees are fitted under codon models either way.
 * --codeml_models (M0, M1a, M2a, M3, M7, M8: src/phylobot/codeml.py) runs PAML codeml on
 * each alignment's best codon combination by BIC (--codeml_trees best, the default) or
 * on every codon combination (all): CODEML per model and starting omega, CODEML_EXPORT
 * per combination (<combination>/codeml/).
 *
 * --branch_site runs only the branch-site test of positive selection (codeml Model A
 * against its null; src/phylobot/branch_site.py) on one codon alignment (--alignment,
 * --seqtype codon) and one tree (--tree) whose foreground branches carry codeml's marks:
 * #1 after a node marks the branch leading to it, $1 every branch of its clade. Results
 * go to <outdir>/branch_site/. The web app starts it from the tree viewer.
 *
 * --properties_only measures property posteriors on ancestors a finished
 * job already published: --properties_dir holds its <aligner>/<tree>/<model>/ files and
 * --property_targets names the combination and node of each ancestor to measure. No
 * reconstruction runs; results go to <outdir>/<aligner>/<tree>/<model>/properties/. The
 * web app starts it from the ancestor viewer, for properties added since the job ran.
 *
 * --uncertainty_only decomposes the variation of a finished job's ancestors into the
 * posterior's own and that of each choice (aligner, tree, model, gap placement):
 * --uncertainty_dir holds the job's published files, and
 * --uncertainty_settings the run's weights and gap versions. Results go to
 * <outdir>/uncertainty/. The web app starts it from the job's Uncertainty page.
 *
 * --simulate spec.json runs only a simulation (src/phylobot/simulate):
 * for each of the spec's seeds, SIM_PREPARE (Python) plans AliSim's runs, SIM_ALISIM
 * runs them in the IQ-TREE container, and SIM_COLLECT writes the truth bundle (input,
 * true alignment, ancestors and tree, provenance) to <outdir>/replicates/<seed>/. A spec
 * on a given tree or root sequence takes it as --sim_tree / --sim_root.
 *
 * --esr runs only extant-sequence reconstruction (src/phylobot/esr.py) on one
 * protein alignment, tree and model, staged in --esr_dir with settings.json (the levels,
 * e0 and/or e1, and the sequences to hold out): ESR_FIT fits the model once, ESR_TIP
 * reconstructs each held-out sequence as if it were an ancestor, and ESR_COLLECT scores
 * them, to <outdir>/esr/. The web app starts it from the ESR page.
 *
 * --species_map / --species_tree reconcile every combination's tree with the species
 * tree (SPECIES_TREE, RECONCILE, RECONCILE_EXPORT). With them,
 * --tree_builders generax adds GeneRax's reconciled tree as a tree source (amino-acid
 * models only): TREE_GENERAX searches from the alignment's IQ-TREE tree (else RAxML-NG's,
 * else the first provided tree), FIT_GENERAX fits it, and TOPOLOGY_TEST runs one AU test
 * per alignment x model over every tree source (<aligner>/topology_tests/<model>.tsv).
 *
 * --stop_after alignments|trees|ancestors (default ancestors) ends the job after that
 * stage: an alignment-only job publishes <outdir>/<aligner>/alignment.fasta (as every
 * job does); a trees-only job each combination's tree, fit and supports (EXPORT_TREE).
 *
 * Data flows as labelled tuples, [meta, files...]. `meta` holds short, stable keys
 * only (aligner, tree), never user text: it names output paths and is part of each
 * task's cache key. Every alignment x tree source x model combination runs
 * independently and is published to <outdir>/<aligner>/<tree>/<model>/.
 *
 * Sources: aligners (mafft, ...) and provided alignments (provided1, ...); tree
 * builders (iqtree, ...) and provided trees (provided1, ...). A provided tree is
 * fitted to each alignment (FIT_IQTREE) rather than searched; --tree_branch_lengths keep
 * (per tree) holds its branch lengths as given. Models: 'modelfinder'
 * (ModelFinder chooses, once per alignment in MODELFINDER, for every tree source) and/or
 * each fixed catalog model (src/phylobot/catalog.py).
 *
 * Pairing rule: never combine separate channels by position (Nextflow pairs by
 * arrival order). Pair by key: `join(failOnMismatch: true, failOnDuplicate: true)`
 * for one-to-one, `combine(by: 0)` for one-to-many.
 */

include { VALIDATE         } from './modules/validate'
include { ALIGN_MAFFT      } from './modules/align_mafft'
include { ALIGN_MUSCLE     } from './modules/align_muscle'
include { ALIGN_CLUSTALO   } from './modules/align_clustalo'
include { ALIGN_PRANK      } from './modules/align_prank'
include { ALIGN_PRANK_CODON } from './modules/align_prank_codon'
include { CODON_ALIGNMENT  } from './modules/codon_alignment'
include { isMixture        } from './modules/mixtures'
include { MODELFINDER      } from './modules/modelfinder_iqtree'
include { MODELFINDER as MODELFINDER_CODON } from './modules/modelfinder_iqtree'
include { TREE_IQTREE      } from './modules/tree_iqtree'
include { TREE_RAXML       } from './modules/tree_raxml'
include { FIT_IQTREE       } from './modules/fit_iqtree'
include { FIT_IQTREE as FIT_CODON_AA } from './modules/fit_iqtree'
include { FIT_IQTREE as FIT_GENERAX } from './modules/fit_iqtree'
include { ANCESTRAL_IQTREE } from './modules/ancestral_iqtree'
include { ROOT_TREE        } from './modules/root_tree'
include { GAPS_ML          } from './modules/gaps_ml'
include { GAPS_PRANK       } from './modules/gaps_prank'
include { GAPS_GRASP       } from './modules/gaps_grasp'
include { EXPORT           } from './modules/export'
include { PROPERTIES       } from './modules/properties'
include { CODEML           } from './modules/codeml'
include { CODEML_EXPORT    } from './modules/codeml_export'
include { BRANCH_SITE_EXPORT } from './modules/branch_site_export'
include { EXPORT_ALIGNMENT } from './modules/export_alignment'
include { EXPORT_TREE      } from './modules/export_tree'
include { SUMMARY          } from './modules/summary'
include { ANCESTOR_MAP     } from './modules/ancestor_map'
include { SPECIES_TREE     } from './modules/species_tree'
include { RECONCILE        } from './modules/reconcile'
include { RECONCILE_EXPORT } from './modules/reconcile_export'
include { TREE_GENERAX     } from './modules/tree_generax'
include { TOPOLOGY_TEST    } from './modules/topology_test'
include { UNCERTAINTY      } from './modules/uncertainty'
include { SIM_PREPARE      } from './modules/sim_prepare'
include { SIM_ALISIM       } from './modules/sim_alisim'
include { SIM_COLLECT      } from './modules/sim_collect'
include { ESR_FIT          } from './modules/esr_fit'
include { ESR_TIP          } from './modules/esr_tip'
include { ESR_COLLECT      } from './modules/esr_collect'

// A parameter given as a comma-separated string (command line) or a list (params file).
def asList(value) {
    if (value == null || value == '' || value == true) {
        return []
    }
    def items = value instanceof List ? value : value.toString().tokenize(',')
    return items.collect { item -> item.toString().trim() }.findAll { item -> item }
}

def parseTools(value, String option, List available) {
    def items = asList(value)
    if (items == ['none']) {
        return []
    }
    if (!items) {
        error "--${option} is empty (available: ${available.join(', ')}, or none)"
    }
    def unknown = items.findAll { item -> !(item in available) }
    if (unknown) {
        error "Unsupported --${option} ${unknown.join(', ')} (available: ${available.join(', ')}, or none)"
    }
    return items.unique()
}

// Catalog model names only (VALIDATE checks them against the catalog); the pattern keeps
// anything else out of paths and commands before validation runs.
def parseModels(boolean codon) {
    def models = asList(params.models)
    // ModelFinder's candidates; unset means --models (as before --candidates existed).
    def candidates = params.candidates == null ? models.findAll { !isMixture(it) } : asList(params.candidates)
    def codonModels = asList(params.codon_models)
    def codonCandidates = asList(params.codon_candidates)
    def codonModelFinder = params.codon_modelfinder ? true : false
    def bad = (models + candidates).findAll { model -> !(model ==~ /[A-Za-z][A-Za-z0-9.]*(\+[A-Z][0-9]{0,2})*/) }
    bad += (codonModels + codonCandidates).findAll { model -> !(model ==~ /[A-Za-z][A-Za-z0-9]*(\+[A-Z][0-9A-Z]{0,3})*/) }
    if (bad) {
        error "Invalid --models, --candidates, --codon_models or --codon_candidates ${bad.unique().join(', ')}"
    }
    if ((codonModels || codonModelFinder) && !codon) {
        error "--codon_models and --codon_modelfinder need coding DNA (--seqtype codon)"
    }
    if (!params.modelfinder && !params.fixed_models && !codonModels && !codonModelFinder) {
        error "Enable --modelfinder and/or --fixed_models, or --codon_modelfinder and/or --codon_models"
    }
    if (params.fixed_models && !models) {
        error "--fixed_models needs --models"
    }
    return [
        models     : models,
        candidates : candidates,
        codonModels: codonModels.unique(),
        codonCandidates: codonCandidates,
        // codon combinations' model keys: codon ModelFinder, then each codon model
        codonKeys  : (codonModelFinder ? ['modelfinder_codon'] : []) + codonModels.unique(),
        // amino-acid combinations' model keys (codon ones: codonModels)
        keys       : (params.modelfinder ? ['modelfinder'] : []) + (params.fixed_models ? models : []),
    ]
}

// A ModelFinder combination's model key (protein or codon ModelFinder).
def isModelFinder(key) {
    return key in ['modelfinder', 'modelfinder_codon']
}

// Which ModelFinder pick a combination uses: its alignment's, of its data type.
def pickKey(meta) {
    return [meta.aligner, meta.codon ? true : false]
}

// The BIC of an IQ-TREE report (the tree step's final fit).
def bic(report) {
    def line = report.readLines().find { text -> text.startsWith('Bayesian information criterion (BIC) score:') }
    if (!line) {
        error "No BIC in ${report}"
    }
    return line.tokenize(' ')[-1] as Double
}

// Branch support replicates: 0 (off) or 1000-10000 (IQ-TREE refuses fewer than 1000).
def checkReplicates(value, String option) {
    def text = value.toString()
    if (!(text ==~ /[0-9]+/) || !((text as Integer) == 0 || (text as Integer) in 1000..10000)) {
        error "--${option} must be 0 (off) or 1000 to 10000 replicates"
    }
}

// `keysValue` (alignments only, --alignment_keys): per file an aligner's key, for an
// alignment inherited from the job a job extends (it keeps the key: VALIDATE checks it),
// or 'provided' for the user's own. Only the user's own count towards max_provided.
def providedFiles(value, labelsValue, String option, keysValue = null) {
    def paths = asList(value)
    def labels = asList(labelsValue)
    def keys = asList(keysValue)
    if (labels && labels.size() != paths.size()) {
        error "--${option}_labels has ${labels.size()} labels for ${paths.size()} files"
    }
    if (keys && keys.size() != paths.size()) {
        error "--${option}_keys has ${keys.size()} keys for ${paths.size()} files"
    }
    def unknown = keys.findAll { key -> !(key in ['provided', 'mafft', 'muscle', 'clustalo', 'prank', 'prank_codon']) }
    if (unknown) {
        error "Unsupported --${option}_keys ${unknown.join(', ')}"
    }
    def own = keys ? keys.count { key -> key == 'provided' } : paths.size()
    if (own > params.max_provided) {
        error "At most ${params.max_provided} files allowed for --${option}"
    }
    return paths.withIndex().collect { path, i ->
        def key = keys && keys[i] != 'provided' ? keys[i] : null
        [file: file(path, checkIfExists: true), label: labels ? labels[i] : null, key: key]
    }
}

// Gap placement methods (src/phylobot/gaps.py), in one place: the inline run takes all
// of them, while PROPERTIES_ONLY keeps only the masked versions, the ones that carry
// IQ-TREE's posteriors and so can be drawn from.
def parseGapMethods(value) {
    return parseTools(value, 'gap_methods', ['parsimony', 'fitch', 'ml', 'prank', 'grasp'])
}

// The gap versions a property posterior can be measured on: PRANK's and GRASP's
// reconstructions carry no posteriors, so nothing can be drawn from them.
def maskedVersions() {
    return ['parsimony', 'fitch', 'ml']
}

// The property keys that may be measured. They reach a command line, so they are checked
// here as well; phylobot.properties.REGISTRY is authoritative and phylobot-properties
// checks them a third time. Both the inline run and PROPERTIES_ONLY come through here.
def parsePropertyKeys(value) {
    if (!value) {
        return []
    }
    return parseTools(value, 'properties', ['charge', 'moment', 'hydropathy_regions', 'disorder', 'esm_likelihood'])
}

// The properties only the phylobot-ml image can compute (torch and model weights;
// phylobot.properties.REGISTRY's `image`). A PROPERTIES task asked for any of them runs
// in that image, a superset of phylobot-tools, so one task still writes one summary.
def mlPropertyKeys() {
    return ['disorder', 'esm_likelihood']
}

workflow {
    if (params.simulate) {
        SIMULATE()
    } else if (params.branch_site) {
        BRANCH_SITE()
    } else if (params.properties_only) {
        PROPERTIES_ONLY()
    } else if (params.uncertainty_only) {
        UNCERTAINTY_ONLY()
    } else if (params.esr) {
        ESR()
    } else {
        PHYLOBOT()
    }
}

// The branch-site test alone: VALIDATE (the codon alignment, the marked tree, the
// control files), CODEML per model and start, BRANCH_SITE_EXPORT.
workflow BRANCH_SITE {
    if (params.seqtype != 'codon') {
        error "--branch_site needs coding DNA (--seqtype codon)"
    }
    if (!(params.genetic_code.toString() ==~ /[0-9]{1,2}/)) {
        error "--genetic_code must be an NCBI genetic code number (1 = standard)"
    }
    def providedAlignments = providedFiles(params.alignment, params.alignment_labels, 'alignment')
    def providedTrees = providedFiles(params.tree, params.tree_labels, 'tree')
    if (params.input || providedAlignments.size() != 1 || providedTrees.size() != 1) {
        error "--branch_site takes one --alignment (a codon alignment) and one --tree with #1 marks, and no --input"
    }
    def manifest = groovy.json.JsonOutput.toJson([
        sequences    : false,
        alignments   : providedAlignments.collect { p -> [filename: p.file.name, label: p.label] },
        trees        : providedTrees.collect { p -> [filename: p.file.name, label: p.label] },
        max_sequences: params.max_sequences,
        max_length   : params.max_length,
        max_provided : params.max_provided,
        seqtype      : 'codon',
        genetic_code : params.genetic_code.toString() as Integer,
        branch_site  : true,
    ])
    VALIDATE(
        channel.value(manifest).collectFile(name: 'manifest.json'),
        channel.value([]),
        channel.value(providedAlignments.collect { p -> p.file }),
        channel.value(providedTrees.collect { p -> p.file }),
        channel.value([]),
        channel.value([]),
    )
    def meta = [aligner: 'provided1', tree: 'provided1', model: 'branch_site']
    def inputs = VALIDATE.out.codon_alignments.flatten()
        .map { alignment -> [meta, alignment] }
        .join(VALIDATE.out.trees.flatten().map { tree -> [meta, tree] }, failOnMismatch: true, failOnDuplicate: true)
    def controls = VALIDATE.out.codeml_controls
        .flatten()
        .map { ctl -> [ctl.parent.name, ctl.baseName, ctl] }
    CODEML(inputs.combine(controls))
    BRANCH_SITE_EXPORT(
        CODEML.out.run
            .groupTuple()
            .join(CODEML.out.control.groupTuple(), failOnMismatch: true, failOnDuplicate: true)
            .join(inputs, failOnMismatch: true, failOnDuplicate: true),
        VALIDATE.out.names.first(),
    )

    VALIDATE.out.versions.mix(CODEML.out.versions)
        .map { file -> file.text }
        .unique()
        .collectFile(name: 'versions.tsv', storeDir: params.outdir, sort: true)
    VALIDATE.out.report
        .collectFile(name: 'validation.json', storeDir: params.outdir)
    VALIDATE.out.sources
        .collectFile(name: 'sources.tsv', storeDir: params.outdir)
}

// Measuring properties on ancestors another job already reconstructed. No new process and
// no new measurement: the existing PROPERTIES runs over staged combination directories
// (--properties_dir), each with an explicit list of the nodes to measure
// (--property_targets, a TSV of combination and node).
//
// The node matching that produced that list is the web tier's (phylobot.analogs needs
// rooted trees and a job-wide view); the pipeline is handed the answer and never repeats
// it. Nothing is validated upstream here, so the staged names are checked against a safe
// pattern before they reach a path or a command line.
workflow PROPERTIES_ONLY {
    def keys = parsePropertyKeys(params.properties)
    if (!keys) {
        error "--properties_only needs --properties"
    }
    if ((params.property_draws as Integer) < 1) {
        error "--property_draws must be at least 1"
    }
    if (!params.properties_dir || !params.property_targets) {
        error "--properties_only needs --properties_dir and --property_targets"
    }
    if (params.input || params.alignment || params.tree) {
        error "--properties_only measures published ancestors, so it takes no --input, --alignment or --tree"
    }
    // Which gap versions to measure: whichever of these the staged files hold
    // (phylobot.pipeline.properties.versions_present skips the rest). PRANK's and
    // GRASP's reconstructions carry no posteriors and cannot be drawn from at all.
    def versions = parseGapMethods(params.gap_methods ?: 'none').findAll { m -> m in maskedVersions() }
    def root = file(params.properties_dir, checkIfExists: true)
    def targets = channel.fromPath(params.property_targets, checkIfExists: true)
        .splitCsv(sep: '\t', header: true)
        .map { row ->
            def combination = (row.combination ?: '').trim()
            def node = (row.node ?: '').trim()
            if (!(combination ==~ /[A-Za-z0-9_.+-]+\/[A-Za-z0-9_.+-]+\/[A-Za-z0-9_.+-]+/)) {
                error "property target: '${combination}' is not an <aligner>/<tree>/<model> combination"
            }
            if (!(node ==~ /[A-Za-z0-9_.-]+/)) {
                error "property target: '${node}' is not an ancestral node name"
            }
            [combination, node]
        }
        .groupTuple()
        .map { combination, nodes ->
            def parts = combination.tokenize('/')
            // type: 'any': files() matches plain files only by default, so anything staged
            // as a directory would be left behind without a word. (Staged today: the
            // version's ancestors and the posterior rows of the nodes asked for.)
            def published = files("${root}/${combination}/*", type: 'any')
            if (!published) {
                error "property target: nothing staged for ${combination}"
            }
            [[aligner: parts[0], tree: parts[1], model: parts[2]], published, nodes.unique().sort().join(',')]
        }
    PROPERTIES(targets, keys.join(','), versions.join(','), mlPropertyKeys())
}

// Where a finished job's ancestral uncertainty comes from: one
// UNCERTAINTY task over its staged results. Node matching, the column frame and the
// weights are all phylobot-uncertainty's, which checks every combination key the staged
// summary names before building a path from it.
workflow UNCERTAINTY_ONLY {
    if (!params.uncertainty_dir || !params.uncertainty_settings) {
        error "--uncertainty_only needs --uncertainty_dir and --uncertainty_settings"
    }
    if (params.input || params.alignment || params.tree) {
        error "--uncertainty_only reads published ancestors, so it takes no --input, --alignment or --tree"
    }
    UNCERTAINTY(
        file(params.uncertainty_dir, checkIfExists: true, type: 'dir'),
        file(params.uncertainty_settings, checkIfExists: true),
    )
}

// Extant-sequence reconstruction alone: ESR_FIT once, ESR_TIP per held-out sequence (by
// the ID ESR_FIT gave it, checked here before it reaches a command line), ESR_COLLECT.
workflow ESR {
    if (!params.esr_dir) {
        error "--esr needs --esr_dir (alignment.fasta, tree.nwk, model.txt, settings.json)"
    }
    if (params.input || params.alignment || params.tree) {
        error "--esr reads --esr_dir, so it takes no --input, --alignment or --tree"
    }
    ESR_FIT(file(params.esr_dir, checkIfExists: true, type: 'dir'))
    def tips = ESR_FIT.out.tips
        .splitText()
        .map { line -> line.trim() }
        .filter { tip -> tip }
        .map { tip ->
            if (!(tip ==~ /t[0-9]{1,6}/)) {
                error "ESR: '${tip}' is not a tip ID"
            }
            tip
        }
    ESR_TIP(ESR_FIT.out.fit.combine(tips))
    ESR_COLLECT(ESR_FIT.out.fit, ESR_TIP.out.result.collect())
}

// A simulation alone: one SIM_PREPARE -> SIM_ALISIM -> SIM_COLLECT chain per seed.
workflow SIMULATE {
    if (params.input || params.alignment || params.tree) {
        error "--simulate takes a spec (and --sim_tree / --sim_root), not --input, --alignment or --tree"
    }
    def spec = file(params.simulate, checkIfExists: true)
    if (spec.extension != 'json') {
        error "--simulate takes a JSON spec (src/phylobot/simulate/spec.py)"
    }
    def data = new groovy.json.JsonSlurper().parse(spec)
    def seeds = data?.simulation?.seeds
    if (!(seeds instanceof List) || !seeds) {
        error "--simulate: the spec has no simulation.seeds"
    }
    def onFile = data.tree?.random == null
    if (onFile != (params.sim_tree as boolean)) {
        error onFile ? "--simulate: this spec needs --sim_tree" : "--simulate: a random-tree spec takes no --sim_tree"
    }
    def fileRoot = data.root?.length == null
    if (fileRoot != (params.sim_root as boolean)) {
        error fileRoot ? "--simulate: this spec needs --sim_root" : "--simulate: a random-root spec takes no --sim_root"
    }
    SIM_PREPARE(
        channel.fromList(seeds),
        channel.value(spec),
        channel.value(params.sim_tree ? file(params.sim_tree, checkIfExists: true) : []),
        channel.value(params.sim_root ? file(params.sim_root, checkIfExists: true) : []),
    )
    SIM_ALISIM(SIM_PREPARE.out.prepared)
    SIM_COLLECT(SIM_ALISIM.out.simulated, channel.value(spec))
}

// The whole pipeline: alignments -> trees -> ancestors (and codeml).
workflow PHYLOBOT {
    // Tools with a module wired in below. Adding a tool: add its key here and a branch.
    def aligners = parseTools(params.aligners, 'aligners', ['mafft', 'muscle', 'clustalo', 'prank', 'prank_codon'])
    def seqtype = params.seqtype ?: 'protein'
    if (!(seqtype in ['protein', 'codon'])) {
        error "--seqtype must be protein or codon (coding DNA)"
    }
    def codon = seqtype == 'codon'
    // An NCBI genetic code number; VALIDATE checks it is one PhyloBot offers.
    if (!(params.genetic_code.toString() ==~ /[0-9]{1,2}/)) {
        error "--genetic_code must be an NCBI genetic code number (1 = standard)"
    }
    if ('prank_codon' in aligners && !codon) {
        error "The prank_codon aligner needs coding DNA (--seqtype codon)"
    }
    def treeBuilders = parseTools(params.tree_builders, 'tree_builders', ['iqtree', 'raxml', 'generax'])
    // Codon combinations' trees from IQ-TREE and RAxML-NG: a codon search and/or fits on
    // the amino-acid trees (the best per tree source, or all of them).
    def codonTrees = parseTools(params.codon_trees ?: 'search', 'codon_trees', ['search', 'amino'])
    def codonAminoTrees = params.codon_amino_trees ?: 'best'
    if (!(codonAminoTrees in ['best', 'all'])) {
        error "--codon_amino_trees must be best or all"
    }
    def providedAlignments = providedFiles(params.alignment, params.alignment_labels, 'alignment', params.alignment_keys)
    def providedTrees = providedFiles(params.tree, params.tree_labels, 'tree')
    // Per provided tree: 'fit' (IQ-TREE re-estimates its branch lengths; the default) or
    // 'keep' (the lengths as given, e.g. a simulation's true ones: -blfix in FIT_IQTREE
    // and ANCESTRAL_IQTREE, modules/branch_lengths.nf; VALIDATE checks every branch has one).
    def treeLengths = asList(params.tree_branch_lengths)
    if (treeLengths && treeLengths.size() != providedTrees.size()) {
        error "--tree_branch_lengths has ${treeLengths.size()} values for ${providedTrees.size()} trees"
    }
    if (treeLengths.any { value -> !(value in ['fit', 'keep']) }) {
        error "--tree_branch_lengths takes fit or keep for each provided tree"
    }
    // Inherited alignments (a job extending another): never also made here.
    def inherited = providedAlignments.findResults { p -> p.key }
    if (inherited.intersect(aligners)) {
        error "--aligners ${inherited.intersect(aligners).join(', ')} already given by --alignment_keys"
    }
    // An extension that adds models or tree sources to the alignments it inherits
    // (src/phylobot/extend.py): the combinations the extended job has already are not
    // run again (--skip_combinations, aligner/tree/model keys), and its ModelFinder picks
    // are reused (--modelfinder_reused, aligner keys; each pick's model.txt and
    // modelfinder.tsv staged in --modelfinder_picks), never recomputed.
    def skipped = asList(params.skip_combinations) as Set
    def badSkips = skipped.findAll { key -> !(key ==~ /[a-z0-9_]+\/[a-z0-9_]+\/[A-Za-z0-9.+_]+/) }
    if (badSkips) {
        error "Invalid --skip_combinations ${badSkips.join(', ')}"
    }
    def reused = asList(params.modelfinder_reused)
    def reusedPicks = reused.collect { aligner ->
        if (!(aligner ==~ /[a-z0-9_]+/) || !params.modelfinder_picks) {
            error "--modelfinder_reused ${aligner} needs its files in --modelfinder_picks"
        }
        def model = file("${params.modelfinder_picks}/${aligner}.model.txt", checkIfExists: true)
        // The pick goes into other steps' command lines (quoted), as MODELFINDER's does.
        if (!(model.text.trim() ==~ /[A-Za-z0-9.+]+/)) {
            error "--modelfinder_picks: ${aligner}'s model does not look like a model"
        }
        [[aligner, false], model, file("${params.modelfinder_picks}/${aligner}.modelfinder.tsv", checkIfExists: true)]
    }
    def stopAfter = params.stop_after ?: 'ancestors'
    if (!(stopAfter in ['alignments', 'trees', 'ancestors'])) {
        error "--stop_after must be alignments, trees or ancestors"
    }
    def withTrees = stopAfter != 'alignments'
    def withAncestors = stopAfter == 'ancestors'
    if (!withTrees && providedTrees) {
        error "Provided trees need --stop_after trees or ancestors"
    }
    def models = withTrees ? parseModels(codon) : [models: [], candidates: [], codonModels: [], codonCandidates: [], codonKeys: [], keys: []]
    checkReplicates(params.bootstrap, 'bootstrap')
    checkReplicates(params.alrt, 'alrt')
    // The outgroup (user sequence names; checked in VALIDATE, passed via the manifest).
    def outgroup = asList(params.outgroup)
    // Gap placement in ancestors (src/phylobot/gaps.py). Unset: fitch with an outgroup,
    // parsimony without (gaps.default_methods); an empty list means none. Only for jobs
    // that reconstruct ancestors.
    def gapDefault = outgroup ? 'fitch' : 'parsimony'
    def gapMethods = !withAncestors ? [] : parseGapMethods(params.gap_methods == null ? gapDefault : (params.gap_methods ?: 'none'))
    // Tier-0 property posteriors. Unset or empty means none, so a
    // hand-run pipeline and the golden scenarios are untouched. The keys are checked
    // here because they reach a command line; phylobot.properties.REGISTRY is
    // authoritative and phylobot-properties checks them again.
    def propertyKeys = !withAncestors ? [] : parsePropertyKeys(params.properties)
    if (propertyKeys && (params.property_draws as Integer) < 1) {
        error "--property_draws must be at least 1"
    }
    def codemlModels = withAncestors ? asList(params.codeml_models) : []
    def codemlTrees = params.codeml_trees ?: 'best'
    if (!(codemlTrees in ['best', 'all'])) {
        error "--codeml_trees must be best or all"
    }
    if (codemlModels && !(models.codonKeys)) {
        error "--codeml_models needs codon models (--codon_modelfinder or --codon_models)"
    }
    if (models.codonKeys && 'amino' in codonTrees) {
        if (!models.keys) {
            error "--codon_trees amino needs amino-acid models (--modelfinder and/or --fixed_models)"
        }
        if (!treeBuilders) {
            error "--codon_trees amino needs --tree_builders iqtree and/or raxml"
        }
    }
    // Reconciliation with a species tree: both files, or neither.
    if ((params.species_map ? true : false) != (params.species_tree ? true : false)) {
        error "Reconciliation needs both --species_map and --species_tree"
    }
    def reconcile = params.species_map ? true : false
    if (reconcile && !withTrees) {
        error "Reconciliation needs trees (--stop_after trees or ancestors)"
    }
    // GeneRax's reconciled tree (a tree source): needs the species tree, and another
    // source's tree to start from (IQ-TREE's, else RAxML-NG's, else the first provided).
    def generax = 'generax' in treeBuilders
    if (generax && !reconcile) {
        error "--tree_builders generax needs a species tree (--species_map and --species_tree)"
    }
    if (generax && !('iqtree' in treeBuilders) && !('raxml' in treeBuilders) && !providedTrees) {
        error "--tree_builders generax needs another tree source to start from (iqtree, raxml or --tree)"
    }
    def generaxStart = 'iqtree' in treeBuilders ? 'iqtree' : 'raxml' in treeBuilders ? 'raxml' : 'provided1'
    if (!(params.polytomies in ['guided', 'arbitrary'])) {
        error "--polytomies must be guided or arbitrary"
    }
    def needRoot = gapMethods.findAll { method -> method in ['prank', 'grasp'] }
    if (needRoot && !outgroup) {
        error "Gap methods ${needRoot.join(', ')} need an --outgroup (their indel models depend on the root)"
    }

    if (!params.input && !providedAlignments) {
        error "Give --input (a protein or, with --seqtype codon, coding-DNA FASTA file) and/or --alignment"
    }
    if (!aligners && !providedAlignments) {
        error "No alignments: --aligners is none and no --alignment was given"
    }
    if (withTrees && !treeBuilders && !providedTrees) {
        error "No trees: --tree_builders is none and no --tree was given"
    }

    // User-controlled names and labels travel only inside this file, never in a command.
    def manifest = groovy.json.JsonOutput.toJson([
        sequences    : params.input ? true : false,
        alignments   : providedAlignments.collect { p -> [filename: p.file.name, label: p.label, key: p.key] },
        trees        : providedTrees.withIndex().collect { p, i ->
            [filename: p.file.name, label: p.label] + (treeLengths && treeLengths[i] == 'keep' ? [keep_lengths: true] : [:])
        },
        max_sequences: params.max_sequences,
        max_length   : params.max_length,
        max_provided : params.max_provided,
        models       : models.models ?: null,
        candidates   : models.candidates ?: null,
        codon_models : models.codonModels ?: null,
        codon_candidates : models.codonCandidates ?: null,
        codon_modelfinder: 'modelfinder_codon' in models.codonKeys,
        codeml_models: codemlModels ?: null,
        raxml        : withTrees && 'raxml' in treeBuilders,
        outgroup     : outgroup,
        seqtype      : seqtype,
        genetic_code : params.genetic_code.toString() as Integer,
        species      : reconcile,
    ])
    VALIDATE(
        channel.value(manifest).collectFile(name: 'manifest.json'),
        params.input ? channel.value(file(params.input, checkIfExists: true)) : channel.value([]),
        channel.value(providedAlignments.collect { p -> p.file }),
        channel.value(providedTrees.collect { p -> p.file }),
        reconcile ? channel.value(file(params.species_map, checkIfExists: true)) : channel.value([]),
        reconcile ? channel.value(file(params.species_tree, checkIfExists: true)) : channel.value([]),
    )
    def names = VALIDATE.out.names.first()  // value channel: reused by every EXPORT
    def sequences = VALIDATE.out.fasta.map { fasta -> [[:], fasta] }  // [meta, fasta]
    def versions = VALIDATE.out.versions
    def modelfinderArgs = VALIDATE.out.modelfinder_args.first()  // value channel
    def codonModelfinderArgs = VALIDATE.out.codon_modelfinder_args.first()
    def modelKeys = channel.fromList(models.keys)

    // Alignments: [meta(aligner), alignment] (protein; for coding DNA, of the translations)
    def provided = VALIDATE.out.alignments
        .flatten()
        .map { alignment -> [[aligner: alignment.baseName], alignment] }
    def aligned = channel.empty()  // from the protein aligners
    if ('mafft' in aligners) {
        ALIGN_MAFFT(sequences)
        aligned = aligned.mix(ALIGN_MAFFT.out.alignment)
        versions = versions.mix(ALIGN_MAFFT.out.versions)
    }
    if ('muscle' in aligners) {
        ALIGN_MUSCLE(sequences)
        aligned = aligned.mix(ALIGN_MUSCLE.out.alignment)
        versions = versions.mix(ALIGN_MUSCLE.out.versions)
    }
    if ('clustalo' in aligners) {
        ALIGN_CLUSTALO(sequences)
        aligned = aligned.mix(ALIGN_CLUSTALO.out.alignment)
        versions = versions.mix(ALIGN_CLUSTALO.out.versions)
    }
    if ('prank' in aligners) {
        ALIGN_PRANK(sequences)
        aligned = aligned.mix(ALIGN_PRANK.out.alignment)
        versions = versions.mix(ALIGN_PRANK.out.versions)
    }

    def alignments = provided.mix(aligned)
    // Codon alignments, [meta(aligner), codon alignment] (coding DNA only), keyed like
    // `alignments`: each aligner's in both forms (CODON_ALIGNMENT), provided ones as given.
    def codonAlignments = channel.empty()
    if (codon) {
        def cds = VALIDATE.out.cds.first()  // value channel
        def codonForm = channel.empty()
        if ('prank_codon' in aligners) {
            ALIGN_PRANK_CODON(VALIDATE.out.cds.map { fasta -> [[:], fasta] })
            codonForm = ALIGN_PRANK_CODON.out.alignment
            versions = versions.mix(ALIGN_PRANK_CODON.out.versions)
        }
        CODON_ALIGNMENT(
            aligned.map { meta, alignment -> [meta, alignment, 'protein'] }
                .mix(codonForm.map { meta, alignment -> [meta, alignment, 'codon'] }),
            cds,
        )
        alignments = provided.mix(CODON_ALIGNMENT.out.protein)
        codonAlignments = VALIDATE.out.codon_alignments
            .flatten()
            .map { alignment -> [[aligner: alignment.baseName], alignment] }
            .mix(CODON_ALIGNMENT.out.codon)
    }

    // Every alignment, published with its original names (the whole result of an
    // alignment-only job); for coding DNA with its codon form.
    def exportedAlignments = codon
        ? alignments.join(codonAlignments, failOnMismatch: true, failOnDuplicate: true)
        : alignments.map { meta, alignment -> [meta, alignment, []] }
    EXPORT_ALIGNMENT(exportedAlignments, names)

    def fits = channel.empty()  // fit.tsv per combination, for SUMMARY
    def exported = channel.empty()  // [meta, published files] per combination, for ANCESTOR_MAP
    if (withTrees) {
        // Trees: [meta(aligner, model, tree), tree, model.txt]; every alignment x every model.
        def alignmentsByModel = alignments
            .combine(modelKeys)
            .map { meta, alignment, model -> [meta + [model: model], alignment] }
        // Codon models: every codon alignment x every codon model key, marked meta.codon.
        if (models.codonKeys) {
            alignmentsByModel = alignmentsByModel.mix(codonAlignments
                .combine(channel.fromList(models.codonKeys))
                .map { meta, alignment, model -> [meta + [model: model, codon: true], alignment] })
        }

        // ModelFinder runs once per alignment (MODELFINDER). Its pick is the model of that
        // alignment's 'modelfinder' combinations for every tree source, so trees from
        // different sources are compared under one model. picks: [aligner, model.txt, report]
        // picks: [[aligner, codon?], model.txt, report]; codon ModelFinder runs on the codon
        // alignments, with its own candidates, for the 'modelfinder_codon' combinations.
        def picks = channel.fromList(reusedPicks)
        // Tree sources, for which of an alignment's ModelFinder combinations still run.
        def sourceKeys = treeBuilders.findAll { t -> t in ['iqtree', 'raxml'] } +
            providedTrees.withIndex().collect { _tree, i -> "provided${i + 1}".toString() }
        def needsPick = { aligner -> sourceKeys.any { t -> !("${aligner}/${t}/modelfinder".toString() in skipped) } }
        if ('modelfinder' in models.keys) {
            MODELFINDER(alignments.filter { meta, _alignment -> !(meta.aligner in reused) && needsPick.call(meta.aligner) }, modelfinderArgs)
            picks = picks.mix(MODELFINDER.out.model.map { meta, model, report -> [[meta.aligner, false], model, report] })
            versions = versions.mix(MODELFINDER.out.versions)
        }
        if ('modelfinder_codon' in models.codonKeys) {
            MODELFINDER_CODON(codonAlignments.map { meta, alignment -> [meta + [codon: true], alignment] }, codonModelfinderArgs)
            picks = picks.mix(MODELFINDER_CODON.out.model.map { meta, model, report -> [[meta.aligner, true], model, report] })
            versions = versions.mix(MODELFINDER_CODON.out.versions)
        }
        // [meta, alignment, picked model]: the pick for a ModelFinder key, [] for a fixed model.
        def modelfinderCombos = alignmentsByModel
            .filter { meta, _alignment -> isModelFinder(meta.model) }
            .map { meta, alignment -> [pickKey(meta), meta, alignment] }
            .combine(picks.map { key, model, _report -> [key, model] }, by: 0)
            .map { _key, meta, alignment, model -> [meta, alignment, model] }
        def fixedCombos = alignmentsByModel
            .filter { meta, _alignment -> !isModelFinder(meta.model) }
            .map { meta, alignment -> [meta, alignment, []] }
        def treeInputs = modelfinderCombos.mix(fixedCombos)
        // Whether a combination still runs here (an extension skips the extended job's).
        def runs = { meta, tree -> !("${meta.aligner}/${tree}/${meta.model}".toString() in skipped) }

        def trees = channel.empty()
        if ('iqtree' in treeBuilders) {
            // Codon models search their own tree only with --codon_trees search.
            TREE_IQTREE(treeInputs.filter { meta, _alignment, _model -> (!meta.codon || 'search' in codonTrees) && runs.call(meta, 'iqtree') })
            trees = trees.mix(TREE_IQTREE.out.tree)
            versions = versions.mix(TREE_IQTREE.out.versions)
        }
        // Topologies that IQ-TREE fits rather than searches (FIT_IQTREE):
        // [meta + tree key, alignment, topology, picked model].
        def fitInputs = channel.empty()
        if ('raxml' in treeBuilders) {
            // RAxML-NG searches; its topology then gets the same model fit as every source.
            def raxmlInputs = treeInputs.filter { meta, _alignment, _model -> !meta.codon && runs.call(meta, 'raxml') }
            TREE_RAXML(raxmlInputs)
            versions = versions.mix(TREE_RAXML.out.versions)
            fitInputs = fitInputs.mix(TREE_RAXML.out.tree
                .join(raxmlInputs, failOnMismatch: true, failOnDuplicate: true)
                .map { meta, tree, alignment, model -> [meta + [tree: 'raxml'], alignment, tree, model] })
        }
        if (providedTrees) {
            // Every alignment x model x provided tree.
            fitInputs = fitInputs.mix(treeInputs
                .combine(VALIDATE.out.trees.flatten())
                .map { meta, alignment, model, tree -> [meta + [tree: tree.baseName], alignment, tree, model] }
                .filter { meta, _alignment, _tree, _model -> runs.call(meta, meta.tree) })
        }
        if ('raxml' in treeBuilders || providedTrees) {
            FIT_IQTREE(fitInputs)
            trees = trees.mix(FIT_IQTREE.out.tree)
            versions = versions.mix(FIT_IQTREE.out.versions)
        }
        if (models.codonKeys && 'amino' in codonTrees) {
            // Amino-acid trees from IQ-TREE and RAxML-NG: [aligner, tree key, tree]. 'best':
            // per alignment and tree source, the amino-acid model of lowest BIC (ties: the
            // first model key), once all its models' trees are in; 'all': every one.
            def aminoTrees = trees
                .filter { meta, _tree, _model, _report -> !meta.codon && meta.tree in ['iqtree', 'raxml'] }
            // One amino-acid model: its tree is the best (and groupTuple(size: 1) never
            // emitted here, Nextflow 26.04.6, so it is not used for that case).
            def bases = codonAminoTrees == 'all'
                ? aminoTrees.map { meta, tree, _model, _report -> [meta.aligner, "${meta.tree}_aa-${meta.model}".toString(), tree] }
                : models.keys.size() == 1
                ? aminoTrees.map { meta, tree, _model, _report -> [meta.aligner, "${meta.tree}_aa".toString(), tree] }
                : aminoTrees
                    .map { meta, tree, _model, report -> [[meta.aligner, meta.tree], [bic(report), meta.model, tree]] }
                    .groupTuple(size: models.keys.size())
                    .map { key, scored ->
                        def best = scored.min { a, b -> a[0] <=> b[0] ?: a[1] <=> b[1] }
                        [key[0], "${key[1]}_aa".toString(), best[2]]
                    }
            // Every codon combination x every amino-acid tree of its alignment.
            def codonFits = treeInputs
                .filter { meta, _alignment, _model -> meta.codon }
                .map { meta, alignment, model -> [meta.aligner, meta, alignment, model] }
                .combine(bases, by: 0)
                .map { _aligner, meta, alignment, model, key, tree -> [meta + [tree: key], alignment, tree, model] }
            FIT_CODON_AA(codonFits)
            trees = trees.mix(FIT_CODON_AA.out.tree)
            versions = versions.mix(FIT_CODON_AA.out.versions)
        }

        // Reconciliation: once every tree is in, the species tree's
        // polytomies are resolved from them (SPECIES_TREE). GeneRax's reconciled tree, a
        // tree source, comes after it and is not among the trees that guide it.
        if (reconcile) {
            SPECIES_TREE(
                VALIDATE.out.species_map,
                VALIDATE.out.species_tree,
                VALIDATE.out.species_names,
                VALIDATE.out.outgroup,
                trees.map { _meta, tree, _model, _report -> tree }.collect(),
            )
        }
        if (generax) {
            // Per alignment x amino-acid model: GeneRax searches from that combination's
            // starting tree (TREE_GENERAX), then IQ-TREE fits the combination's model on
            // the topology (FIT_GENERAX), as for RAxML-NG.
            // [meta(aligner, model), alignment, picked model, starting tree, its model.txt]
            def generaxInputs = treeInputs
                .filter { meta, _alignment, _model -> !meta.codon }
                .join(trees
                    .filter { meta, _tree, _model, _report -> !meta.codon && meta.tree == generaxStart }
                    .map { meta, tree, model, _report -> [meta.findAll { k, _v -> k != 'tree' }, tree, model] },
                    failOnMismatch: true, failOnDuplicate: true)
            TREE_GENERAX(
                generaxInputs.map { meta, alignment, _picked, tree, model -> [meta, alignment, tree, model] },
                SPECIES_TREE.out.tree.first(),
                SPECIES_TREE.out.mapping.first(),
            )
            versions = versions.mix(TREE_GENERAX.out.versions)
            FIT_GENERAX(TREE_GENERAX.out.tree
                .join(generaxInputs, failOnMismatch: true, failOnDuplicate: true)
                .map { meta, tree, _generaxModel, alignment, picked, _start, _model -> [meta + [tree: 'generax'], alignment, tree, picked] })
            trees = trees.mix(FIT_GENERAX.out.tree)

            // Every tree source's topology under each alignment x amino-acid model, in one
            // AU test (TOPOLOGY_TEST): is the reconciled tree acceptable to the sequences?
            def sources = treeBuilders.size() + providedTrees.size()
            def candidates = trees
                .filter { meta, _tree, _model, _report -> !meta.codon }
                .map { meta, tree, model, _report -> [meta.findAll { k, _v -> k != 'tree' }, [meta.tree, tree, model]] }
                .groupTuple(size: sources)
                .map { meta, entries ->
                    def sorted = entries.sort { a, b -> a[0] <=> b[0] }
                    [meta, sorted[0][2], sorted.collect { e -> e[0] }, sorted.collect { e -> e[1] }]
                }
            TOPOLOGY_TEST(treeInputs
                .filter { meta, _alignment, _model -> !meta.codon }
                .map { meta, alignment, _picked -> [meta, alignment] }
                .join(candidates, failOnMismatch: true, failOnDuplicate: true))
            versions = versions.mix(TOPOLOGY_TEST.out.versions)
        }

        // IQ-TREE reports, keyed by the full meta, for EXPORT: [meta, the tree step's report
        // (fit statistics), ModelFinder's report (candidate ranking; [] for a fixed model)].
        def rankings = trees
            .filter { meta, _tree, _model, _report -> isModelFinder(meta.model) }
            .map { meta, _tree, _model, _report -> [pickKey(meta), meta] }
            .combine(picks.map { key, _model, report -> [key, report] }, by: 0)
            .map { _key, meta, report -> [meta, report] }
            .mix(trees
                .filter { meta, _tree, _model, _report -> !isModelFinder(meta.model) }
                .map { meta, _tree, _model, _report -> [meta, []] })
        def reports = trees
            .map { meta, _tree, _model, report -> [meta, report] }
            .join(rankings, failOnMismatch: true, failOnDuplicate: true)

        // Each tree rejoins the alignment it was built from (one alignment, many trees):
        // the codon alignment for a codon model. `proteinInputs` has the protein one for
        // every combination (gap placement and exports work per protein column).
        def alignmentsByType = alignments
            .map { meta, alignment -> [[meta.aligner, false], alignment] }
            .mix(codonAlignments.map { meta, alignment -> [[meta.aligner, true], alignment] })
        def asrInputs = trees
            .map { meta, tree, model, _report -> [[meta.aligner, meta.codon ? true : false], meta, tree, model] }
            .combine(alignmentsByType, by: 0)
            .map { _key, meta, tree, model, alignment -> [meta, alignment, tree, model] }
        def proteinInputs = trees
            .map { meta, tree, model, _report -> [meta.aligner, meta, tree, model] }
            .combine(alignments.map { meta, alignment -> [meta.aligner, alignment] }, by: 0)
            .map { _key, meta, tree, model, alignment -> [meta, alignment, tree, model] }
        // [meta, codon alignment] for codon combinations, [meta, []] otherwise (EXPORT).
        def codonExport = asrInputs
            .map { meta, alignment, _tree, _model -> [meta, meta.codon ? alignment : []] }
        // Each combination's tree with the node names the job publishes, for reconciliation:
        // the reconstruction's (Node1...), or the search tree in a trees-only job (named
        // by the export as EXPORT_TREE names it).
        def publishedTrees = proteinInputs.map { meta, _alignment, tree, _model -> [meta, tree] }

        if (withAncestors) {
            ANCESTRAL_IQTREE(asrInputs)
            versions = versions.mix(ANCESTRAL_IQTREE.out.versions)
            publishedTrees = ANCESTRAL_IQTREE.out.ancestors.map { meta, asrTree, _state -> [meta, asrTree] }

            // Gap placement steps that need another program run per combination (ml, prank,
            // grasp), each on that combination's alignment and tree. Each yields
            // [meta, fileA, fileB]; a method not chosen yields empty placeholders.
            def perCombination = ANCESTRAL_IQTREE.out.ancestors
                .join(proteinInputs, failOnMismatch: true, failOnDuplicate: true)
                .map { meta, asrTree, _state, alignment, _tree, model -> [meta, alignment, asrTree, model] }
            def placeholder = asrInputs.map { meta, _alignment, _tree, _model -> [meta, [], []] }

            // With an outgroup, each combination's tree is rooted on it: published, used by Fitch
            // (in EXPORT), and given to PRANK and GRASP, whose indel models depend on the root.
            def rooted = placeholder
            if (outgroup) {
                ROOT_TREE(
                    ANCESTRAL_IQTREE.out.ancestors.map { meta, asrTree, _state -> [meta, asrTree] },
                    VALIDATE.out.outgroup.first(),
                )
                rooted = ROOT_TREE.out.rooted
            }
            def withRooted = perCombination
                .join(rooted, failOnMismatch: true, failOnDuplicate: true)
                .map { meta, alignment, _asrTree, model, rootedTree, _report -> [meta, alignment, rootedTree, model] }

            def mlGaps = placeholder
            def prankGaps = placeholder
            def graspGaps = placeholder
            if ('ml' in gapMethods) {
                GAPS_ML(perCombination.map { meta, alignment, tree, _model -> [meta, alignment, tree] })
                mlGaps = GAPS_ML.out.gaps
            }
            if ('prank' in gapMethods) {
                GAPS_PRANK(withRooted.map { meta, alignment, tree, _model -> [meta, alignment, tree] })
                prankGaps = GAPS_PRANK.out.ancestors
            }
            if ('grasp' in gapMethods) {
                GAPS_GRASP(withRooted)
                graspGaps = GAPS_GRASP.out.ancestors
            }

            // One-to-one by the full meta: each reconstruction with its own inputs.
            def exportInputs = ANCESTRAL_IQTREE.out.ancestors
                .join(proteinInputs, failOnMismatch: true, failOnDuplicate: true)
                .join(reports, failOnMismatch: true, failOnDuplicate: true)
                .join(mlGaps, failOnMismatch: true, failOnDuplicate: true)
                .join(prankGaps, failOnMismatch: true, failOnDuplicate: true)
                .join(graspGaps, failOnMismatch: true, failOnDuplicate: true)
                .join(rooted, failOnMismatch: true, failOnDuplicate: true)
                .join(codonExport, failOnMismatch: true, failOnDuplicate: true)
                .join(ANCESTRAL_IQTREE.out.sites, failOnMismatch: true, failOnDuplicate: true)
            EXPORT(exportInputs, names, gapMethods.join(','))

            // Tier-0 property posteriors: pure Python on each
            // combination's published ancestors, so it needs only EXPORT's output.
            if (propertyKeys) {
                PROPERTIES(
                    EXPORT.out.results.map { meta, published -> [meta, published, ''] },
                    propertyKeys.join(','),
                    gapMethods.join(','),
                    mlPropertyKeys(),
                )
            }

            fits = EXPORT.out.fit
            exported = EXPORT.out.results

            // PAML codeml on the codon combinations: each alignment's best by BIC (from
            // the exported fits) or every one, x each model's starting omegas.
            if (codemlModels) {
                def codonFits = EXPORT.out.fit.filter { meta, _fit -> meta.codon }
                def chosen = codemlTrees == 'all'
                    ? codonFits.map { meta, _fit -> [meta, true] }
                    : codonFits
                        .map { meta, fit -> [meta.aligner, meta, fit.readLines()[1].tokenize('\t')[5] as Double] }
                        .groupTuple()
                        .map { _aligner, metas, bics -> [metas[bics.indexOf(bics.min())], true] }
                def codemlInputs = chosen
                    .join(ANCESTRAL_IQTREE.out.ancestors, failOnDuplicate: true)
                    .join(asrInputs, failOnDuplicate: true)
                    .map { meta, _chosen, asrTree, _state, alignment, _tree, _model -> [meta, alignment, asrTree] }
                def controls = VALIDATE.out.codeml_controls
                    .flatten()
                    .map { ctl -> [ctl.parent.name, ctl.baseName, ctl] }
                CODEML(codemlInputs
                    .combine(controls)
                    .map { meta, alignment, tree, model, start, ctl -> [meta, alignment, tree, model, start, ctl] })
                versions = versions.mix(CODEML.out.versions)
                CODEML_EXPORT(CODEML.out.run
                    .groupTuple()
                    .join(codemlInputs.map { meta, _alignment, tree -> [meta, tree] }, failOnMismatch: true, failOnDuplicate: true),
                    names)
            }
        } else {
            EXPORT_TREE(
                proteinInputs.join(reports, failOnMismatch: true, failOnDuplicate: true),
                names,
                VALIDATE.out.outgroup.first(),
            )
            fits = EXPORT_TREE.out.fit
            exported = EXPORT_TREE.out.results
        }

        // GeneRax reconciles each combination's tree with the species tree (RECONCILE) and
        // the events go onto the published tree (RECONCILE_EXPORT).
        if (reconcile) {
            RECONCILE(proteinInputs, SPECIES_TREE.out.tree.first(), SPECIES_TREE.out.mapping.first())
            versions = versions.mix(RECONCILE.out.versions)
            RECONCILE_EXPORT(
                RECONCILE.out.reconciliation.join(publishedTrees, failOnMismatch: true, failOnDuplicate: true),
                names,
                SPECIES_TREE.out.nodes.first(),
                !withAncestors,
            )
        }
    }

    // Inherited alignments first: they were the extended job's aligners.
    def primaryAligner = (inherited + aligners) ? (inherited + aligners)[0] : 'provided1'
    def primaryTree = treeBuilders ? treeBuilders[0] : 'provided1'
    def primary = withTrees ? "${primaryAligner}/${primaryTree}/${(models.keys + models.codonKeys)[0]}" : ''
    // One row per combination: its keys, then its fit.tsv values (header skipped). An
    // alignment-only job has no combinations: an empty table.
    def fitRows = withTrees
        ? fits
            .map { meta, fit -> "${meta.aligner}\t${meta.tree}\t${meta.model}\t" + fit.readLines()[1] }
            .collectFile(name: 'fits.tsv', newLine: true, sort: true)
        : channel.of('').collectFile(name: 'fits.tsv', newLine: true)
    def alignmentKeys = alignments
        .map { meta, _alignment -> meta.aligner }
        .collect()
        .map { keys -> keys.sort().join(',') }
    SUMMARY(fitRows, primary, VALIDATE.out.sources, alignmentKeys, stopAfter)

    // Every combination's ancestors by content-derived IDs (phylobot.analogs): one line
    // per combination with its key, whether its tree is outgroup-rooted, and the published
    // tree. Nothing is exported without trees, so then it does not run.
    def mapTrees = exported
        .map { meta, files ->
            def published = [files].flatten()
            def rooted = published.find { f -> f.name == 'tree.rooted.nwk' }
            def tree = rooted ?: published.find { f -> f.name == 'tree.nwk' }
            "${meta.aligner}/${meta.tree}/${meta.model}\t${rooted ? 'rooted' : 'unrooted'}\t${tree.text.trim()}"
        }
        .collectFile(name: 'trees.tsv', newLine: true, sort: true)
    ANCESTOR_MAP(mapTrees)

    versions
        .map { file -> file.text }
        .unique()
        .collectFile(name: 'versions.tsv', storeDir: params.outdir, sort: true)
    VALIDATE.out.report
        .collectFile(name: 'validation.json', storeDir: params.outdir)
    VALIDATE.out.sources
        .collectFile(name: 'sources.tsv', storeDir: params.outdir)
}
