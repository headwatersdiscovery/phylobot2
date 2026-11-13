// Site-profile mixtures (catalog.MIXTURE_VARIANTS; a test pins this rule to it). IQ-TREE keeps
// C20's published weights unless told to fit them, so every step that fits a mixture passes
// -mwopt. Other models' commands are unchanged. Only fixed models can be mixtures: ModelFinder
// never picks one, so a meta.model of 'modelfinder' is never one.
def isMixture(model) {
    return model ==~ /.*\+C[0-9]+(\+.*)?/
}

def mixtureFlag(meta) {
    return isMixture(meta.model) ? ' -mwopt' : ''
}
