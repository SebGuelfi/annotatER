setwd("/home/sguelfi/projects/R/annotatER/")
library(devtools)
load_all()

## set the path and take the names of the samples
## to create the proper path in ubuntu I have to do mount the neuroscience path like this.
## sshfs skgtmgu@wise.cs.ucl.ac.uk:/SAN/neuroscience /home/sguelfi/neuroscience/

path <- "/home/sguelfi/neuroscience/ukbec/analysis/hipp/2ndPassSTAR/"
list.samples = list.dirs(path = path,recursive = F)
names(list.samples) <- basename(list.samples)
list.samples[1:length(list.samples)] <- paste0(list.samples,"/align_star.bamSJ.out.tab")
## list.samples = the list of samples to collect
## minSamples = filter split reads that do not have this min of number of samples

tmp.table <- NULL
minSamples = 5  ## in percentage

STARSplitRead <- loadSplitReads(list.samples = list.samples,minSamples = 5)

## Separate counts and create the junction ID so that can be mapped back
junctions <- list()
junctions[["annotation"]] <- cbind(junID=1:nrow(STARSplitRead),STARSplitRead[,c(1:4,7)])
junctions[["counts"]] <- cbind(junID=1:nrow(STARSplitRead),STARSplitRead[,8:ncol(STARSplitRead)])

#load("~/projects/R/hipp/data/expression/splitReads/splitReads.filtered.rda")

GTFPath <- "/data/references/GTF/Homo_sapiens.GRCh38.87.gtf"

## annotate the junctions
junctions[["annotation"]] <- annotateSplitReads(GTFPath,junctions[["annotation"]])







