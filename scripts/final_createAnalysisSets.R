f.get.simple.analysis.sets <- function() {
  out <- list(
    both_vs_none = list(
      modelFunction = "lm",
      distributionFamily = "none",
      targetVars = c("MSR", "MBD", "MER", "MEV"),
      formulaParts = c("~ seqDepth + Pot + isLupin + correct_An + SymboGroup + Species"),
      alternative_FandP = list("Pot" = "Species", "correct_An" = "Species", "SymboGroup" = "Species", "__Symbogroup" = "Species"),
      combineContrasts = list("__Symbogroup" = c("correct_An", "SymboGroup"))
    ),
    AMF_vs_noAMF = list(
      modelFunction = "lm",
      distributionFamily = "none",
      targetVars = c("MSR", "MBD", "MER", "MEV"),
      formulaParts = c("~ seqDepth + Pot + isLupin + correct_AR + SymboGroup + Species"),
      alternative_FandP = list("Pot" = "Species", "correct_AR" = "Species", "SymboGroup" = "Species", "__Symbogroup" = "Species"),
      combineContrasts = list("__Symbogroup" = c("correct_AR", "SymboGroup"))
    ),
    Rhizo_vs_noRhizo = list(
      modelFunction = "lm",
      distributionFamily = "none",
      targetVars = c("MSR", "MBD", "MER", "MEV"),
      formulaParts = c("~ seqDepth + Pot + isLupin + correct_nn + SymboGroup + Species"),
      alternative_FandP = list("Pot" = "Species", "correct_nn" = "Species", "SymboGroup" = "Species", "__Symbogroup" = "Species"),
      combineContrasts = list("__Symbogroup" = c("correct_nn", "SymboGroup"))
    ),
    mainFactor = list(
      modelFunction = "lm",
      distributionFamily = "none",
      targetVars = c("MSR", "MBD", "MER", "MEV"),
      formulaParts = c("~ seqDepth + Pot + isLupin + SymboGroup + Species"),
      alternative_FandP = list("Pot" = "Species", "SymboGroup" = "Species"),
      combineContrasts = list()
    )
  )
  return(out)
}

f.get.species.contrasts <- function(terms) {
  cat("WARNING: THE TERMS HAVE TO USE UP THE FACTOR SPECIES!\n")
  out <- list(
    noBlockWithPotAgainstSpecies = list(
      modelFunction = "lm",
      distributionFamily = "none",
      targetVars = c("MSR", "MBD", "MER", "MEV"),
      formulaParts = c("~ seqDepth + Pot + SymboGroup + ", paste0(terms, collapse = " + ")),
      alternative_FandP = list("Pot" = "__Species", "SymboGroup" = "__Species"),
      combineContrasts = list("__Species" = c(terms))
    ),
    noBlockNoPotAgainstSpecies = list(
      modelFunction = "lm",
      distributionFamily = "none",
      targetVars = c("MSR", "MBD", "MER", "MEV"),
      formulaParts = c("~ seqDepth + SymboGroup + ", paste0(terms, collapse = " + ")),
      alternative_FandP = list("SymboGroup" = "__Species"),
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
    SymboGroup = list(
      modelFunction = "lm",
      distributionFamily = "none",
      targetVars = c("irrelevant"),
      formulaParts = c("~ Pot + SymboGroup"),
      alternative_FandP = list(),
      combineContrasts = list()
    ),
    SymboGroupWithBlock = list(
      modelFunction = "lm",
      distributionFamily = "none",
      targetVars = c("irrelevant"),
      formulaParts = c("~ Block + Pot + SymboGroup"),
      alternative_FandP = list(),
      combineContrasts = list()
    )#,
    #within_AMF_nonRhizo = list(
    #  modelFunction = "lm",
    #  distributionFamily = "none",
    #  targetVars = c("irrelevant"),
    #  formulaParts = c("~ Block + Pot + SymboGroup + within_AMF_nonRhizo"), # interaction must be removed, it's otherwise not full rank
    #  alternative_FandP = list(),
    #  combineContrasts = list()
    #),
    #within_AMF_Rhizo = list(
    #  modelFunction = "lm",
    #  distributionFamily = "none",
    #  targetVars = c("irrelevant"),
    #  formulaParts = c("~ Block + Pot + SymboGroup + within_AMF_Rhizo"), # interaction must be removed, it's otherwise not full rank
    #  alternative_FandP = list(),
    #  combineContrasts = list()
    #),
    #within_nonAMF_nonRhizo = list(
    #  modelFunction = "lm",
    #  distributionFamily = "none",
    #  targetVars = c("irrelevant"),
    #  formulaParts = c("~ Block + Pot + SymboGroup + within_nonAMF_nonRhizo"), # interaction must be removed, it's otherwise not full rank
    #  alternative_FandP = list(),
    #  combineContrasts = list()
    #)
  )
  return(out)
}

f.get.simple.analysis.sets.for.da.high.vs.low <- function() { # removed block
  out <- list(
    Rhizobia = list(
      modelFunction = "lm",
      distributionFamily = "none",
      targetVars = c("irrelevant"),
      formulaParts = c("~ Pot + RhizobiaBinary"),
      alternative_FandP = list(),
      combineContrasts = list()
    ),
    AMF = list(
      modelFunction = "lm",
      distributionFamily = "none",
      targetVars = c("irrelevant"),
      formulaParts = c("~ Pot + AMFbinary"),
      alternative_FandP = list(),
      combineContrasts = list()
    ),
    RhizobiaWithinAMF = list(
      modelFunction = "lm",
      distributionFamily = "none",
      targetVars = c("irrelevant"),
      formulaParts = c("~ Pot + RhizobiaWithinAMF"),
      alternative_FandP = list(),
      combineContrasts = list()
    ),
    AMFwithinNonRhizo = list(
      modelFunction = "lm",
      distributionFamily = "none",
      targetVars = c("irrelevant"),
      formulaParts = c("~ Pot + AMFwithinNonRhizo"),
      alternative_FandP = list(),
      combineContrasts = list()
    ),
    bothVsNone = list(
      modelFunction = "lm",
      distributionFamily = "none",
      targetVars = c("irrelevant"),
      formulaParts = c("~ Pot + bothVsNone"),
      alternative_FandP = list(),
      combineContrasts = list()
    )
  )
  return(out)
}

f.get.species.contrasts.da  <- function(x) { # removed block
  if (grepl("Lupin", x, ignore.case = TRUE)) {
    f.print.message("WARNING: Lupin cannot be fit with SymboGroup. Skipping SymboGroup!")
    out <- list(
      modelFunction = "lm",
      distributionFamily = "none",
      targetVars = c("irrelevant"),
      formulaParts = c("~ Pot + ", x),
      alternative_FandP = list(),
      combineContrasts = list()
    )
  } else {
    out <- list(
      modelFunction = "lm",
      distributionFamily = "none",
      targetVars = c("irrelevant"),
      formulaParts = c("~ Pot + SymboGroup + ", x),
      alternative_FandP = list(),
      combineContrasts = list()
    )
  }
  return(out)
}



