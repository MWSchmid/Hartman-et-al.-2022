f.get.simple.analysis.sets <- function() {
  out <- list(
    withBlockAndPotAgainstSpecies = list(
      modelFunction = "lm",
      distributionFamily = "none",
      targetVars = c("MSR", "MBD", "MER", "MEV"),
      formulaParts = c("~ Block + Pot + Rhizobia + AMF + Rhizobia:AMF + Species"),
      alternative_FandP = list("Pot" = "Species", "Rhizobia" = "Species", "AMF" = "Species", "Rhizobia:AMF" = "Species"),
      combineContrasts = list()
    ),
    withBlockNoPotAgainstSpecies = list(
      modelFunction = "lm",
      distributionFamily = "none",
      targetVars = c("MSR", "MBD", "MER", "MEV"),
      formulaParts = c("~ Block + Rhizobia + AMF + Rhizobia:AMF + Species"),
      alternative_FandP = list("Rhizobia" = "Species", "AMF" = "Species", "Rhizobia:AMF" = "Species"),
      combineContrasts = list()
    ),
    noBlockWithPotAgainstSpecies = list(
      modelFunction = "lm",
      distributionFamily = "none",
      targetVars = c("MSR", "MBD", "MER", "MEV"),
      formulaParts = c("~ Pot + Rhizobia + AMF + Rhizobia:AMF + Species"),
      alternative_FandP = list("Pot" = "Species", "Rhizobia" = "Species", "AMF" = "Species", "Rhizobia:AMF" = "Species"),
      combineContrasts = list()
    ),
    noBlockNoPotAgainstSpecies = list(
      modelFunction = "lm",
      distributionFamily = "none",
      targetVars = c("MSR", "MBD", "MER", "MEV"),
      formulaParts = c("~ Rhizobia + AMF + Rhizobia:AMF + Species"),
      alternative_FandP = list("Rhizobia" = "Species", "AMF" = "Species", "Rhizobia:AMF" = "Species"),
      combineContrasts = list()
    ),
    withBlockAndPotAgainstSpeciesFirstAMF = list(
      modelFunction = "lm",
      distributionFamily = "none",
      targetVars = c("MSR", "MBD", "MER", "MEV"),
      formulaParts = c("~ Block + Pot + AMF + Rhizobia + AMF:Rhizobia + Species"),
      alternative_FandP = list("Pot" = "Species", "Rhizobia" = "Species", "AMF" = "Species", "AMF:Rhizobia" = "Species"),
      combineContrasts = list()
    ),
    withBlockNoPotAgainstSpeciesFirstAMF = list(
      modelFunction = "lm",
      distributionFamily = "none",
      targetVars = c("MSR", "MBD", "MER", "MEV"),
      formulaParts = c("~ Block + AMF + Rhizobia + AMF:Rhizobia + Species"),
      alternative_FandP = list("Rhizobia" = "Species", "AMF" = "Species", "AMF:Rhizobia" = "Species"),
      combineContrasts = list()
    ),
    noBlockWithPotAgainstSpeciesFirstAMF = list(
      modelFunction = "lm",
      distributionFamily = "none",
      targetVars = c("MSR", "MBD", "MER", "MEV"),
      formulaParts = c("~ Pot + AMF + Rhizobia + AMF:Rhizobia + Species"),
      alternative_FandP = list("Pot" = "Species", "Rhizobia" = "Species", "AMF" = "Species", "AMF:Rhizobia" = "Species"),
      combineContrasts = list()
    ),
    noBlockNoPotAgainstSpeciesFirstAMF = list(
      modelFunction = "lm",
      distributionFamily = "none",
      targetVars = c("MSR", "MBD", "MER", "MEV"),
      formulaParts = c("~ AMF + Rhizobia + AMF:Rhizobia + Species"),
      alternative_FandP = list("Rhizobia" = "Species", "AMF" = "Species", "AMF:Rhizobia" = "Species"),
      combineContrasts = list()
    )#,
    # Note that species first does not work because Rhizobia and AMF are out
    # withBlockAndPotSpeciesFirst = list(
    #   modelFunction = "lm",
    #   distributionFamily = "none",
    #   targetVars = c("MSR", "MBD", "MER", "MEV"),
    #   formulaParts = c("~ Block + Pot + Species + Rhizobia + AMF"),
    #   alternative_FandP = list(),
    #   combineContrasts = list()
    # ),
    # withBlockNoPotSpeciesFirst = list(
    #   modelFunction = "lm",
    #   distributionFamily = "none",
    #   targetVars = c("MSR", "MBD", "MER", "MEV"),
    #   formulaParts = c("~ Block + Species + Rhizobia + AMF"),
    #   alternative_FandP = list("Rhizobia" = "Species", "AMF" = "Species"),
    #   combineContrasts = list()
    # ),
    # noBlockWithPotSpeciesFirst = list(
    #   modelFunction = "lm",
    #   distributionFamily = "none",
    #   targetVars = c("MSR", "MBD", "MER", "MEV"),
    #   formulaParts = c("~ Pot + Species + Rhizobia + AMF"),
    #   alternative_FandP = list(),
    #   combineContrasts = list()
    # ),
    # noBlockNoPotSpeciesFirst = list(
    #   modelFunction = "lm",
    #   distributionFamily = "none",
    #   targetVars = c("MSR", "MBD", "MER", "MEV"),
    #   formulaParts = c("~ Species + Rhizobia + AMF"),
    #   alternative_FandP = list(),
    #   combineContrasts = list()
    # )
  )
  return(out)
}

f.get.species.contrasts <- function(terms) {
  cat("WARNING: THE TERMS HAVE TO USE UP THE FACTOR SPECIES!\n")
  out <- list(
    withBlockAndPotAgainstSpecies = list(
      modelFunction = "lm",
      distributionFamily = "none",
      targetVars = c("MSR", "MBD", "MER", "MEV"),
      formulaParts = c("~ Block + Pot + Rhizobia + AMF + Rhizobia:AMF + ", paste0(terms, collapse = " + ")),
      alternative_FandP = list("Pot" = "__Species", "Rhizobia" = "__Species", "AMF" = "__Species", "Rhizobia:AMF" = "__Species"),
      combineContrasts = list("__Species" = c(terms))
    ),
    withBlockNoPotAgainstSpecies = list(
      modelFunction = "lm",
      distributionFamily = "none",
      targetVars = c("MSR", "MBD", "MER", "MEV"),
      formulaParts = c("~ Block + Rhizobia + AMF + Rhizobia:AMF + ", paste0(terms, collapse = " + ")),
      alternative_FandP = list("Rhizobia" = "__Species", "AMF" = "__Species", "Rhizobia:AMF" = "__Species"),
      combineContrasts = list("__Species" = c(terms))
    ),
    noBlockWithPotAgainstSpecies = list(
      modelFunction = "lm",
      distributionFamily = "none",
      targetVars = c("MSR", "MBD", "MER", "MEV"),
      formulaParts = c("~ Pot + Rhizobia + AMF + Rhizobia:AMF + ", paste0(terms, collapse = " + ")),
      alternative_FandP = list("Pot" = "__Species", "Rhizobia" = "__Species", "AMF" = "__Species", "Rhizobia:AMF" = "__Species"),
      combineContrasts = list("__Species" = c(terms))
    ),
    noBlockNoPotAgainstSpecies = list(
      modelFunction = "lm",
      distributionFamily = "none",
      targetVars = c("MSR", "MBD", "MER", "MEV"),
      formulaParts = c("~ Rhizobia + AMF + Rhizobia:AMF + ", paste0(terms, collapse = " + ")),
      alternative_FandP = list("Rhizobia" = "__Species", "AMF" = "__Species", "Rhizobia:AMF" = "__Species"),
      combineContrasts = list("__Species" = c(terms))
    ),
    withBlockAndPotAgainstSpeciesFirstAMF = list(
      modelFunction = "lm",
      distributionFamily = "none",
      targetVars = c("MSR", "MBD", "MER", "MEV"),
      formulaParts = c("~ Block + Pot + AMF + Rhizobia + AMF:Rhizobia + ", paste0(terms, collapse = " + ")),
      alternative_FandP = list("Pot" = "__Species", "Rhizobia" = "__Species", "AMF" = "__Species", "AMF:Rhizobia" = "__Species"),
      combineContrasts = list("__Species" = c(terms))
    ),
    withBlockNoPotAgainstSpeciesFirstAMF = list(
      modelFunction = "lm",
      distributionFamily = "none",
      targetVars = c("MSR", "MBD", "MER", "MEV"),
      formulaParts = c("~ Block + AMF + Rhizobia + AMF:Rhizobia + ", paste0(terms, collapse = " + ")),
      alternative_FandP = list("Rhizobia" = "__Species", "AMF" = "__Species", "AMF:Rhizobia" = "__Species"),
      combineContrasts = list("__Species" = c(terms))
    ),
    noBlockWithPotAgainstSpeciesFirstAMF = list(
      modelFunction = "lm",
      distributionFamily = "none",
      targetVars = c("MSR", "MBD", "MER", "MEV"),
      formulaParts = c("~ Pot + AMF + Rhizobia + AMF:Rhizobia + ", paste0(terms, collapse = " + ")),
      alternative_FandP = list("Pot" = "__Species", "Rhizobia" = "__Species", "AMF" = "__Species", "AMF:Rhizobia" = "__Species"),
      combineContrasts = list("__Species" = c(terms))
    ),
    noBlockNoPotAgainstSpeciesFirstAMF = list(
      modelFunction = "lm",
      distributionFamily = "none",
      targetVars = c("MSR", "MBD", "MER", "MEV"),
      formulaParts = c("~ AMF + Rhizobia + AMF:Rhizobia + ", paste0(terms, collapse = " + ")),
      alternative_FandP = list("Rhizobia" = "__Species", "AMF" = "__Species", "AMF:Rhizobia" = "__Species"),
      combineContrasts = list("__Species" = c(terms))
    )
  )
  return(out)
}

f.get.simple.analysis.sets.for.da <- function() {
  out <- list(
    #Rhizobia = list(
    #  modelFunction = "lm",
    #  distributionFamily = "none",
    #  targetVars = c("irrelevant"),
    #  formulaParts = c("~ Block + Pot + Rhizobia"),
    #  alternative_FandP = list(),
    #  combineContrasts = list()
    #),
    #AMF = list(
    #  modelFunction = "lm",
    #  distributionFamily = "none",
    #  targetVars = c("irrelevant"),
    #  formulaParts = c("~ Block + Pot + AMF"),
    #  alternative_FandP = list(),
    #  combineContrasts = list()
    #),
    Block = list(
      modelFunction = "lm",
      distributionFamily = "none",
      targetVars = c("irrelevant"),
      formulaParts = c("~ Block"),
      alternative_FandP = list(),
      combineContrasts = list()
    ),
    Pot = list(
      modelFunction = "lm",
      distributionFamily = "none",
      targetVars = c("irrelevant"),
      formulaParts = c("~ Pot"),
      alternative_FandP = list(),
      combineContrasts = list()
    ),
    PotAfterBlock = list(
      modelFunction = "lm",
      distributionFamily = "none",
      targetVars = c("irrelevant"),
      formulaParts = c("~ Block + Pot"),
      alternative_FandP = list(),
      combineContrasts = list()
    ),
    Rhizobia_X_AMF = list(
      modelFunction = "lm",
      distributionFamily = "none",
      targetVars = c("irrelevant"),
      formulaParts = c("~ Block + Pot + Rhizobia + AMF + Rhizobia:AMF"),
      alternative_FandP = list(),
      combineContrasts = list()
    ),
    within_AMF_nonRhizo = list(
      modelFunction = "lm",
      distributionFamily = "none",
      targetVars = c("irrelevant"),
      formulaParts = c("~ Block + Pot + Rhizobia + AMF + within_AMF_nonRhizo"), # interaction must be removed, it's otherwise not full rank
      alternative_FandP = list(),
      combineContrasts = list()
    ),
    within_AMF_Rhizo = list(
      modelFunction = "lm",
      distributionFamily = "none",
      targetVars = c("irrelevant"),
      formulaParts = c("~ Block + Pot + Rhizobia + AMF + within_AMF_Rhizo"), # interaction must be removed, it's otherwise not full rank
      alternative_FandP = list(),
      combineContrasts = list()
    ),
    within_nonAMF_nonRhizo = list(
      modelFunction = "lm",
      distributionFamily = "none",
      targetVars = c("irrelevant"),
      formulaParts = c("~ Block + Pot + Rhizobia + AMF + within_nonAMF_nonRhizo"), # interaction must be removed, it's otherwise not full rank
      alternative_FandP = list(),
      combineContrasts = list()
    )
  )
  return(out)
}

f.get.simple.analysis.sets.for.da.high.vs.low <- function() {
  out <- list(
    Rhizobia = list(
      modelFunction = "lm",
      distributionFamily = "none",
      targetVars = c("irrelevant"),
      formulaParts = c("~ Block + Pot + RhizobiaBinary"),
      alternative_FandP = list(),
      combineContrasts = list()
    ),
    AMF = list(
      modelFunction = "lm",
      distributionFamily = "none",
      targetVars = c("irrelevant"),
      formulaParts = c("~ Block + Pot + AMFbinary"),
      alternative_FandP = list(),
      combineContrasts = list()
    )
  )
  return(out)
}

f.get.simple.analysis.sets.for.da.high.vs.low.species.last <- function() { # does not work!
  out <- list(
    Rhizobia = list(
      modelFunction = "lm",
      distributionFamily = "none",
      targetVars = c("irrelevant"),
      formulaParts = c("~ Block + Pot + RhizobiaBinary + SpeciesChar"),
      alternative_FandP = list(),
      combineContrasts = list()
    ),
    AMF = list(
      modelFunction = "lm",
      distributionFamily = "none",
      targetVars = c("irrelevant"),
      formulaParts = c("~ Block + Pot + AMFbinary + SpeciesChar"),
      alternative_FandP = list(),
      combineContrasts = list()
    )
  )
  return(out)
}

f.get.species.contrasts.da  <- function(x) {
  out <- list(
      modelFunction = "lm",
      distributionFamily = "none",
      targetVars = c("irrelevant"),
      formulaParts = c("~ Block + Pot  + Rhizobia + AMF + ", x),
      alternative_FandP = list(),
      combineContrasts = list()
  )
  return(out)
}









