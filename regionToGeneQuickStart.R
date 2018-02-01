#' ---
#' title: "Annotation of expressed regions using split reads"
#' author: "Sebastian Guelfi & David Zhang"
#' date: "January 26th, 2018"
#' ---

#' Quick example of how genomic regions can be annotated using the **XXX**

#' First step is to load the junction definitions. In this example we load the
#' transcript definition from GTEX, generated from the **recount** resource.
#' Source: (XXX)

#+ echo=F
library(data.table)
splitReadsInfo <- fread("/data/recount/GTEx/SRP012682.junction_id_with_transcripts.bed")
## We reformat the GTEx to obtain the junID
splitReadsInfo[,"junID" := tstrsplit(V4, "|", fixed=TRUE,keep=1)]
setDF(splitReadsInfo)


#' We load the split reads counts per tissue. The counts for the split reads can
#' be downloaded from here
#' Now, we take as example the junction count expression for putamen
splitReadCounts <- read.csv("/data/recount/GTEx/byTissue/datatable/Brain-Putamen_basalganglia/Brain-Putamen_basalganglia.csv",check.names=FALSE)
#' The count expression data.frame contain n+1 columns where n is the number of
#' samples. The extra column contains the unique ID for the junctions. The name
#' of the unique ID column must be 'junID'

## we select only the OMIM genes
OMIMGenes <- read.csv("/home/dzhang/projects/OMIM/one_time_files_for_seb/OMIM_filtered_uniq_ens_ids_151017.csv",header = F)

#' We set the percentage of the number of samples that have to pass the counts
#' to 5%
system.time(splitReadTable <- filterSplitReads(splitReadCounts,splitReadsInfo, 5))

#' Load the path for the GTF reference file which will be use to annotate the
#' juctions
GTFPath <- "/data/references/GTF/Homo_sapiens.GRCh38.87.gtf"

#' Here we add one to the start and stop position of the UCSC, this is to match
#' the UCSC coordinates to the ENSEMBL GTF reference file, UCSC coordinates
#' start at postion 0, while ENSEMBL coordinates start at position 1. If a UCSC
#' GTF reference file is used this step is no needed.
splitReadTable$start <- splitReadTable$start +1
splitReadTable$stop <- splitReadTable$stop +1

### take a sample of 10000 split reads

#splitReadTable <- head(splitReadTable,10000)
## we select only the first 6 columns that correspond to the coordinates of the
## split reads
splitReadTable <- splitReadTable[,1:6]
head(splitReadTable)

system.time(splitReadTable <- annotateSplitReads(GTFPath,splitReadTable,geneList = OMIMGenes$V1))

### load
load("/home/dzhang/projects/OMIM/results/intron_intergenic_regions_for_Seb_package_test/ERs_overlapping_OMIM_max_gap_1Mb_intron_intergenic_PUTM.rda")

library(SummarizedExperiment)
regions <- rowRanges(ERs_overlapping_OMIM_max_gap_1Mb_intron_intergenic_PUTM)

annotated.regions <- annotatERJunction(regions,splitReadTable)








