#!/usr/bin/env Rscript

metricToUse <- "MICe"
rDir <- "/media/mwschmid/myData/MWSchmid/MarcelVanDerHeijden_rhizospheres/GitIgnore_results_ITS_semiRarefied/OTU/"
subDir <- "coreMicrobiomes"
infileName <- "network_all_values.csv.adj.sig.gz"

##########################################################################################
### arguments from commandline
myarg <- commandArgs()
argPos <- grep("--args", myarg, fixed = TRUE)
metricToUse <- as.character(myarg[argPos+1])
rDir <- as.character(myarg[argPos+2])
subDir <- as.character(myarg[argPos+3])
infileName <- as.character(myarg[argPos+4])

##########################################################################################
### libraries
suppressPackageStartupMessages({
  library("igraph")
  library("poweRlaw")
  library("RColorBrewer")
  library("gplots")
  library("data.table")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapper.R")
  source("/media/mwschmid/myData/MWSchmid/Development/R/RNAseqWrapperAds.R")
})

f.print.message <- function(...) { cat("=== ", format(Sys.time(), "%Y %b %d %X"), paste("===", ...,"\n")) }

pcname <- system('uname -n',intern=T)
numCoresPerMachine <- c("nuke" = 4, "styx" = 15, "marc-IEU" = 7, "piftel" = 30)
numCoresAvailable <- numCoresPerMachine[pcname]

##########################################################################################
### load data
networkDir <- file.path(rDir, subDir)
dir.create(networkDir, showWarnings = FALSE)
infileName <- file.path(rDir, infileName)
tempFile <- file.path("/home/marc/hurzliBurzli.gz")
file.copy(infileName, tempFile)
if (file.exists("/home/marc/hurzliBurzli")) {
  f.print.message("WARNING: DON'T RUN TWO INSTANCES OF THIS SCRIPT IN PARALLEL!!!")
  file.remove("/home/marc/hurzliBurzli")
}
system(paste0("gunzip ", tempFile))
#myData <- read.table(gsub(".gz", "", tempFile), sep = '\t', header = TRUE)
myData <- fread(gsub(".gz", "", tempFile), sep = '\t', data.table = FALSE)
myData <- subset(myData, FDR_TICe < 0.05) # should already be the case
myDataSimple <- myData[,c("varX", metricToUse, "varY", paste0("FDR_", metricToUse))]
colnames(myDataSimple) <- c("otuA", "strength", "otuB", "FDR")
myDataSimple <- subset(myDataSimple, FDR < 0.05) # also strength should be significant

write.table(myDataSimple, file.path(networkDir, paste0(metricToUse, "_forCytoscape.txt")), row.names = FALSE, quote = FALSE, sep = '\t')

write.table(myDataSimple[,c("otuA", "otuB", "strength")], file.path(networkDir, paste0(metricToUse, "_forIGRAPH.txt")), sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)

f.detect.communities.and.roles.in.network(file.path(networkDir, paste0(metricToUse, "_forIGRAPH.txt")), networkDir, paste0(metricToUse, '_'))

# load the communities and roles
communities <- read.csv(file.path(networkDir, paste0(metricToUse, "_communitiesAndRoles.csv")), stringsAsFactors = FALSE, header = TRUE)
rownames(communities) <- communities$geneID; colnames(communities)[1] <- "OTU"
memberCounts <- table(communities$community); memberCounts
write.csv(communities, file.path(networkDir, paste0(metricToUse, "_communitiesAndRoles.csv")), quote = FALSE, row.names = FALSE)

##########################################################################################
### check the network
f.check.undirected.coexpression.network(file.path(networkDir, paste0(metricToUse, "_forIGRAPH.txt")), networkDir, metricToUse, simplifyNetwork = TRUE, checkFit = FALSE, numCores = numCoresAvailable)







