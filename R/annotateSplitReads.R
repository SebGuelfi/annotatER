
## this functions load the split reads from STAR and filter those split reads that are not detected in at least minSamples parameter, also I minNumber of
## of counts can be set.
#' @param path list of junction file samples
#' @param the percentage of samples that have to detect the split reads samples, otherwise split read is filtered
#' @parma minCounts: the minimum number of counts to be detected to count for that split reads
loadSplitReads  <- function(list.samples,minSamples,minCounts=1)
{
    ## collecting the split reads information per sample
    splitReads.list <- list()
    for(j in 1:length(list.samples))
    {
        ## filter based on counts
        tmp.SR <- read.delim(file=list.samples[j],header = F ,sep = "\t",
                             col.names = c("chr", "start", "stop", "strand", "intronMotif", "inAnnotation", "uniquelyMapped","multiMapped", "maxSpliOverhang"))
        tmp.SR <- tmp.SR[(tmp.SR[,"uniquelyMapped"]>= minCounts),]
        ## assign it to the list
        splitReads.list[[names(list.samples)[j]]] <- tmp.SR

        cat("\r ",paste0(round((j/length(list.samples)),digits = 4)*100,"%"), " Samples loaded" )

    }
    cat("\n")
    tmp.table <- NULL
    for(i in 1:length(splitReads.list))
    {
        if(i==1)
        {
            tmp <- splitReads.list[[i]]
            tmp.table <- tmp[,c(1:6)]
            tmp.table$countsSamples <- 1
            rm(tmp)

        }else{

            tmp2 <- splitReads.list[[i]]
            tmp.table$countsSamples[which(as.character(tmp.table$chr)%in%as.character(tmp2$chr)
                                          & as.character(tmp.table$start)%in%as.character(tmp2$start)
                                          & as.character(tmp.table$stop)%in%as.character(tmp2$stop))] <- tmp.table[which(as.character(tmp.table$chr)%in%as.character(tmp2$chr)
                                                                                                                         & as.character(tmp.table$start)%in%as.character(tmp2$start)
                                                                                                                         & as.character(tmp.table$stop)%in%as.character(tmp2$stop)),"countsSamples"] + 1
            if(nrow(tmp2[-which(as.character(tmp2$chr)%in%as.character(tmp.table$chr) & as.character(tmp2$start)%in%as.character(tmp.table$start) & as.character(tmp2$stop)%in%as.character(tmp.table$stop)),c(1:3,6)])>0){
                tmp.table <- rbind(tmp.table,cbind(tmp2[-which(as.character(tmp2$chr)%in%as.character(tmp.table$chr)
                                                               & as.character(tmp2$start)%in%as.character(tmp.table$start)
                                                               & as.character(tmp2$stop)%in%as.character(tmp.table$stop)),c(1:6)],countsSamples=1))
            }
            rm(tmp2)
        }

        cat("\r ",paste0(round((i/length(splitReads.list)),digits = 4)*100,"%"), " Samples processed" )
    }
    rm(i,j)
    cat("\n")
    ## filter based on the number of samples that have the split reads
    tot.tmp <- nrow(tmp.table)
    message(tot.tmp, " total number of split reads detected")
    tmp.table <- tmp.table[which(tmp.table$countsSamples>=(length(list.samples)*0.01*minSamples)),]
    message(paste0(round(100-(nrow(tmp.table)/tot.tmp*100),digits = 2),"%"), " split reads filtered")
    message(paste0(round(((nrow(tmp.table)-sum(tmp.table$inAnnotation))/nrow(tmp.table))*100,digits = 2),"%"), " novel split reads identified")

    cat("\n")
    for(i in names(splitReads.list))
    {
        tmp2 <- splitReads.list[[i]]
        dim(tmp.table)
        tmp.table <- merge(tmp.table,tmp2[,c("chr","start","stop","uniquelyMapped")],by=c("chr","start","stop"),all.x = T)
        colnames(tmp.table)[colnames(tmp.table)=="uniquelyMapped"] <- i
        tmp.table[is.na(tmp.table[,i]),i]<-0
        rm(tmp2)
        cat("\r ",paste0(round((match(i,names(splitReads.list))/length(splitReads.list)),digits = 4)*100,"%"), " Samples processed" )
    }
    cat("\n Completed")
    return(tmp.table)
}


## this function filters the split read
#' @param splitReadCounts The table with the split reads counts per sample
#' @param splitReadInfo the split read table information containing the
#' information about the split reads coordinates
#' @param minSamples the minimum number of samples that have to pass the
#' minCounts threshold to be retain
#' @return minCounts the minimum number of counts that are going to be used as
#' threshold (1 by default)
filterSplitReads <- function(splitReadCounts,splitReadsInfo, minSamples, minCounts=1)
{
    ## filter split reads that do not have at least minCounts split reads counts in at least minSamples number of samples
    idx <- which(rowSums(splitReadCounts[,-which(colnames(splitReadCounts)=="junID")]>=minCounts,na.rm = T)>(ncol(splitReadCounts)-1)*0.01*minSamples)
    splitReads.final <- cbind(splitReadsInfo[match(splitReadCounts[idx,"junID"],splitReadsInfo$junID),c("junID","V1","V2","V3","V6")],
                              cbind(rowSums(splitReadCounts[idx,-which(colnames(splitReadCounts)=="junID")]>=minCounts,na.rm = T),
                                    splitReadCounts[idx,-which(colnames(splitReadCounts)=="junID")]))
    colnames(splitReads.final)[2:6] <- c("chr","start","stop","strand","countsSamples")
    return(splitReads.final)
}


#' @param GTFPath the path to the GTF file that is use to annotate the split
#' reads
#' @param splitReadTable the split read table to be annotated
#' @param format the source of the split reads either GTEX or STAR
#' (GTEX by default)
#' @return dataframe returns a data.frame object with the annotated split reads
annotateSplitReads <- function(GTFPath,splitReadTable,geneList=NULL)
{
    library(data.table)
    library(tidyverse)
    library(refGenome)
    # change to directory where you downloaded GTF file
    #setwd("/home/sguelfi/neuroscience/WT_BRAINEAC/hipp/salmon/salmonReferences/GTF/")
    setwd(dirname(GTFPath))
    # create ensemblGenome object for storing Ensembl genomic annotation data
    ens <- ensemblGenome()
    # read GTF file into ensemblGenome object
    basedir(ens) <- dirname(GTFPath)
    message(paste(Sys.time(),"Loading the GTF..."))
    read.gtf(ens, "Homo_sapiens.GRCh38.87.gtf")
    message(paste(Sys.time(),"GTF Loaded"))

    jens <- getSpliceTable(ens)
    jens <- as.data.frame(jens@ev$gtf)

    if(!is.null(geneList))
    {
        jens <- jens %>% filter(as.character(gene_id)%in%as.character(geneList))
        message(paste(Sys.time(),length(unique(jens$gene_id)),"genes retained"))
        rm(geneList)
    }

    jens$strand <- as.character(jens$strand)


    jens$seqid <- gsub("chr","",as.character(jens$seqid))
    splitReadTable$chr <- gsub("chr","",as.character(splitReadTable$chr))

    # head(jens)
    # id seqid    lstart      lend    rstart      rend         gene_id gene_name strand   transcript_id   lexid   rexid transcript_biotype
    # 1     7 127588345 127588565 127589083 127589163 ENSG00000004059      ARF5      + ENST00000000233 1007522 1007525     protein_coding
    # 2     7 127589083 127589163 127589485 127589594 ENSG00000004059      ARF5      + ENST00000000233 1007525 1007527     protein_coding
    # 3     7 127589485 127589594 127590066 127590137 ENSG00000004059      ARF5      + ENST00000000233 1007527 1007529     protein_coding
    # 4     7 127590066 127590137 127590963 127591088 ENSG00000004059      ARF5      + ENST00000000233 1007529 1007531     protein_coding
    # 5     7 127590963 127591088 127591213 127591705 ENSG00000004059      ARF5      + ENST00000000233 1007531 1007533     protein_coding
    # 6    12   8940365   8941940   8942416   8942542 ENSG00000003056      M6PR      - ENST00000000412 1566634 1566632     protein_coding


    ## -1 the start and +1 the stop of the junctions to match the splice site positions
    ## get the coordinates
    splitReadTable$start <- splitReadTable$start -1
    splitReadTable$stop <- splitReadTable$stop +1

    ## convert the strand from numeric to "+,-,*" because this is what the ensembl reference uses.
    if(any(names(table(splitReadTable$strand))==1) | any(names(table(splitReadTable$strand))==2) | any(names(table(splitReadTable$strand))==0))
    {
        splitReadTable$strand[splitReadTable$strand==1] <- "+"
        splitReadTable$strand[splitReadTable$strand==2] <- "-"
        splitReadTable$strand[splitReadTable$strand==0] <- "*"
    }


    message(paste(Sys.time(), 'Getting acceptor sites annotation'))
    ## get the acceptor for the forward strand
    splitReadTable.forward <- splitReadTable
    splitReadTable.forward.tmp <- left_join(splitReadTable.forward,
                                                             jens[as.character(jens$strand)=="+",c("seqid","lstart","strand","transcript_id")],
                                                             by=c("chr"="seqid","stop"="lstart","strand"="strand"))

    splitReadTable.forward <- left_join(splitReadTable.forward,
                                                             jens[as.character(jens$strand)=="+",c("seqid","rstart","strand","transcript_id")],
                                                             by=c("chr"="seqid","stop"="rstart","strand"="strand"))

    ## get the acceptor for the reverse strand
    splitReadTable.reverse <- splitReadTable
    splitReadTable.reverse.tmp <- left_join(splitReadTable.reverse,
                                            jens[as.character(jens$strand)=="-",c("seqid","lend","strand","transcript_id")],
                                            by=c("chr"="seqid","start"="lend","strand"="strand"))

    splitReadTable.reverse <- left_join(splitReadTable.reverse,
                                            jens[as.character(jens$strand)=="-",c("seqid","rend","strand","transcript_id")],
                                            by=c("chr"="seqid","start"="rend","strand"="strand"))


    ## merge the acceptor for forward and reverse strand
    splitReadTable.acceptor <- rbind(splitReadTable.forward,splitReadTable.forward.tmp,splitReadTable.reverse.tmp,splitReadTable.reverse)
    rm(splitReadTable.forward,splitReadTable.forward.tmp,splitReadTable.reverse.tmp,splitReadTable.reverse)
    ## filter duplicates
    splitReadTable.acceptor <- splitReadTable.acceptor[!duplicated(splitReadTable.acceptor[,c("junID","transcript_id")]),]
    setDT(splitReadTable.acceptor)
    ## collapse the transcripts that have the same location of the split site
    splitReadTable.acceptor <- splitReadTable.acceptor %>% group_by(junID) %>% summarise(transcript_id = paste0(na.omit(transcript_id), collapse=","))
    # head(splitReadTable.acceptor)
    # junId   transcript_id
    # 1              NA
    # 2              NA
    # 3 ENST00000488147
    # 4 ENST00000488147
    # 5 ENST00000488147
    # 6 ENST00000488147

    ## set back to data.frame
    setDF(splitReadTable.acceptor)

    ## assign the transcript to the acceptor column
    splitReadTable[match(splitReadTable.acceptor$junID,splitReadTable$junID),"acceptor"] <- splitReadTable.acceptor$transcript_id
    rm(splitReadTable.acceptor)

    ## get the donor
    message(paste(Sys.time(), 'Getting donor sites annotation'))
    splitReadTable.forward <- splitReadTable

    ## get the donor for the forward strand
    splitReadTable.forward.tmp <- left_join(splitReadTable.forward,
                                            jens[as.character(jens$strand)=="+",c("seqid","lend","strand","transcript_id")],
                                            by=c("chr"="seqid","start"="lend","strand"="strand"))

    splitReadTable.forward <- left_join(splitReadTable.forward,
                                            jens[as.character(jens$strand)=="+",c("seqid","rend","strand","transcript_id")],
                                            by=c("chr"="seqid","start"="rend","strand"="strand"))

    ## get the donor for the reverse strand
    splitReadTable.reverse <- splitReadTable
    splitReadTable.reverse.tmp <- left_join(splitReadTable.reverse,
                                            jens[as.character(jens$strand)=="-",c("seqid","lstart","strand","transcript_id")],
                                            by=c("chr"="seqid","stop"="lstart","strand"="strand"))

    splitReadTable.reverse <- left_join(splitReadTable.reverse,
                                            jens[as.character(jens$strand)=="-",c("seqid","rstart","strand","transcript_id")],
                                            by=c("chr"="seqid","stop"="rstart","strand"="strand"))

    splitReadTable.donor <- rbind(splitReadTable.forward,splitReadTable.forward.tmp,splitReadTable.reverse.tmp,splitReadTable.reverse)
    rm(splitReadTable.forward,splitReadTable.forward.tmp,splitReadTable.reverse.tmp,splitReadTable.reverse)
    ## filter duplicates
    splitReadTable.donor <- splitReadTable.donor[!duplicated(splitReadTable.donor[,c("junID","transcript_id")]),]
    setDT(splitReadTable.donor)
    splitReadTable.donor <- splitReadTable.donor %>% group_by(junID) %>% summarise(transcript_id = paste0(na.omit(transcript_id), collapse=","))
    setDF(splitReadTable.donor)

    splitReadTable[match(splitReadTable.donor$junID,splitReadTable$junID),"donor"] <- splitReadTable.donor$transcript_id
    rm(splitReadTable.donor)
    # head(splitReadTable)
    #   junId chr start  stop strand intronMotif inAnnotation countsSamples        acceptor           donor
    # 1     1   1 14829 14930      -           2            1            12
    # 2     2   1 14829 14970      -           2            1           100
    # 3     3   1 15038 15796      -           2            1            59 ENST00000488147 ENST00000488147
    # 4     4   1 15947 16607      -           2            1            15 ENST00000488147 ENST00000488147
    # 5     5   1 16765 16858      -           2            1            27 ENST00000488147 ENST00000488147
    # 6     6   1 17055 17233      -           2            1            46 ENST00000488147 ENST00000488147

    message(paste(Sys.time(), 'Getting junction sites annotation'))

    splitReadTable.junction <-left_join(splitReadTable,
              jens[,c("seqid","lend","rstart","strand","transcript_id")],
              by=c("chr"="seqid","start"="lend", "stop"="rstart","strand"="strand"))


    ## filter duplicates
    splitReadTable.junction <- splitReadTable.junction[!duplicated(splitReadTable.junction[,c("junID","transcript_id")]),]
    setDT(splitReadTable.junction)
    splitReadTable.junction <- splitReadTable.junction %>% group_by(junID) %>% summarise(transcript_id = paste0(na.omit(transcript_id), collapse=","))
    setDF(splitReadTable.junction)
    splitReadTable[match(splitReadTable.junction$junID,splitReadTable$junID),"junction"] <- splitReadTable.junction$transcript_id
    rm(splitReadTable.junction)
    #head(splitReadTable)

    # junId chr start  stop strand intronMotif inAnnotation countsSamples        acceptor           donor        junction
    # 1     1   1 14829 14930      -           2            1            12
    # 2     2   1 14829 14970      -           2            1           100
    # 3     3   1 15038 15796      -           2            1            59 ENST00000488147 ENST00000488147 ENST00000488147
    # 4     4   1 15947 16607      -           2            1            15 ENST00000488147 ENST00000488147 ENST00000488147
    # 5     5   1 16765 16858      -           2            1            27 ENST00000488147 ENST00000488147 ENST00000488147
    # 6     6   1 17055 17233      -           2            1            46 ENST00000488147 ENST00000488147 ENST00000488147

    #save(splitReadTable,file="~/projects/R/hipp/data/expression/splitReads/jun.table.ann.tmp.rda")
    library(GenomicRanges)
    message(paste(Sys.time(),"Annotating split reads with no exact boundaries"))
    splitReadTable$precBoundDonor <- FALSE
    splitReadTable$precBoundAcceptor <- FALSE
    ## add information about correct boundaries
    splitReadTable[which(splitReadTable$junction !=""),c("precBoundDonor","precBoundAcceptor")] <- TRUE
    splitReadTable[which(splitReadTable$donor!=""),c("precBoundDonor")] <- TRUE
    splitReadTable[which(splitReadTable$acceptor!=""),c("precBoundAcceptor")] <- TRUE

    ## In this section I check whether the split reads lie in coding regions of the transcript, but not on splice site.
    ## forward donor

    splitReadTable.notPreBoun.forward <- splitReadTable[which((!splitReadTable$precBoundDonor) & splitReadTable$strand=="+"),]
    #head(splitReadTable.notPreBoun)

    jens.forward <- jens[jens$strand=="+",]

    juncti.GR <- GRanges(c(as.character(splitReadTable.notPreBoun.forward$chr)),
                         IRanges(start=c(splitReadTable.notPreBoun.forward$start),
                                 end = c(splitReadTable.notPreBoun.forward$start)),
                         strand = c(splitReadTable.notPreBoun.forward$strand),
                         juncId=c(splitReadTable.notPreBoun.forward$junID))

    jens.GR <- GRanges(c(as.character(jens.forward$seqid),as.character(jens.forward$seqid)),
                       IRanges(start=c(jens.forward$lstart,jens.forward$rstart),
                               end = c(jens.forward$lend,jens.forward$rend)),
                       strand = c(as.character(jens.forward$strand),as.character(jens.forward$strand)),
                       juncId=c(jens.forward$id,jens.forward$id))

    tabOver <- as.data.frame(findOverlaps(juncti.GR,jens.GR,ignore.strand=FALSE))
    tabOver <- cbind(tabOver,transcripts=jens[match(mcols(jens.GR[tabOver$subjectHits])[,"juncId"],jens$id),"transcript_id"])
    ## remove duplicates
    tabOver <- tabOver[!duplicated(tabOver[,c("queryHits","transcripts")]),]
    if(nrow(tabOver)>0)
    {
        setDT(tabOver)
        tabOver <- tabOver %>% group_by(queryHits) %>% mutate(transcripts = str_c(na.omit(transcripts), collapse=","))
        setDF(tabOver)
        splitReadTable[match(mcols(juncti.GR[tabOver$queryHits])[,"juncId"],splitReadTable$junID),"donor"] <- tabOver$transcripts
    }

    rm(splitReadTable.notPreBoun.forward,tabOver,jens.GR,juncti.GR)

    ## donor reverse
    splitReadTable.notPreBoun.reverse <- splitReadTable[which((!splitReadTable$precBoundDonor) &splitReadTable$strand=="-"),]
    jens.reverse <- jens[jens$strand=="-",]

    juncti.GR <- GRanges(c(as.character(splitReadTable.notPreBoun.reverse$chr)),
                         IRanges(start=c(splitReadTable.notPreBoun.reverse$stop),
                                 end = c(splitReadTable.notPreBoun.reverse$stop)),
                         strand = c(splitReadTable.notPreBoun.reverse$strand),
                         juncId=c(splitReadTable.notPreBoun.reverse$junID))

    jens.GR <- GRanges(c(as.character(jens.reverse$seqid),as.character(jens.reverse$seqid)),
                       IRanges(start=c(jens.reverse$lstart,jens.reverse$rstart),
                               end = c(jens.reverse$lend,jens.reverse$rend)),
                       strand = c(as.character(jens.reverse$strand),as.character(jens.reverse$strand)),
                       juncId=c(jens.reverse$id,jens.reverse$id))

    tabOver <- as.data.frame(findOverlaps(juncti.GR,jens.GR,ignore.strand=FALSE))
    tabOver <- cbind(tabOver,transcripts=jens[match(mcols(jens.GR[tabOver$subjectHits])[,"juncId"],jens$id),"transcript_id"])
    ## remove duplicates
    tabOver <- tabOver[!duplicated(tabOver[,c("queryHits","transcripts")]),]
    if(nrow(tabOver)>0)
    {
        setDT(tabOver)
        tabOver <- tabOver %>% group_by(queryHits) %>% summarise(transcripts = str_c(na.omit(transcripts), collapse=","))
        setDF(tabOver)
        splitReadTable[match(mcols(juncti.GR[tabOver$queryHits])[,"juncId"],splitReadTable$junID),"donor"] <- tabOver$transcripts
    }
    rm(splitReadTable.notPreBoun.reverse,tabOver,jens.GR,juncti.GR)

    ## Acceptor forward

    splitReadTable.notPreBoun.forward <- splitReadTable[which((!splitReadTable$precBoundAcceptor) & splitReadTable$strand=="+"),]

    juncti.GR <- GRanges(c(as.character(splitReadTable.notPreBoun.forward$chr)),
                         IRanges(start=c(splitReadTable.notPreBoun.forward$stop),
                                 end = c(splitReadTable.notPreBoun.forward$stop)),
                         strand = c(splitReadTable.notPreBoun.forward$strand),
                         juncId=c(splitReadTable.notPreBoun.forward$junID))

    jens.GR <- GRanges(c(as.character(jens.forward$seqid),as.character(jens.forward$seqid)),
                       IRanges(start=c(jens.forward$lstart,jens.forward$rstart),
                               end = c(jens.forward$lend,jens.forward$rend)),
                       strand = c(as.character(jens.forward$strand),as.character(jens.forward$strand)),
                       juncId=c(jens.forward$id,jens.forward$id))

    tabOver <- as.data.frame(findOverlaps(juncti.GR,jens.GR,ignore.strand=FALSE))
    tabOver <- cbind(tabOver,transcripts=jens[match(mcols(jens.GR[tabOver$subjectHits])[,"juncId"],jens$id),"transcript_id"])
    ## remove duplicates

    tabOver <- tabOver[!duplicated(tabOver[,c("queryHits","transcripts")]),]
    if(nrow(tabOver)>0)
    {
        setDT(tabOver)
        tabOver <- tabOver %>% group_by(queryHits) %>% summarise(transcripts = str_c(na.omit(transcripts), collapse=","))
        setDF(tabOver)
        splitReadTable[match(mcols(juncti.GR[tabOver$queryHits])[,"juncId"],splitReadTable$junID),"acceptor"] <- tabOver$transcripts
    }

    rm(splitReadTable.notPreBoun.forward,tabOver,jens.GR,juncti.GR,jens.forward)

    ##Acceptor reverse
    splitReadTable.notPreBoun.reverse <- splitReadTable[which((!splitReadTable$precBoundDonor) &splitReadTable$strand=="-"),]

    juncti.GR <- GRanges(c(as.character(splitReadTable.notPreBoun.reverse$chr)),
                         IRanges(start=c(splitReadTable.notPreBoun.reverse$start),
                                 end = c(splitReadTable.notPreBoun.reverse$start)),
                         strand = c(splitReadTable.notPreBoun.reverse$strand),
                         juncId=c(splitReadTable.notPreBoun.reverse$junID))

    jens.GR <- GRanges(c(as.character(jens.reverse$seqid),as.character(jens.reverse$seqid)),
                       IRanges(start=c(jens.reverse$lstart,jens.reverse$rstart),
                               end = c(jens.reverse$lend,jens.reverse$rend)),
                       strand = c(as.character(jens.reverse$strand),as.character(jens.reverse$strand)),
                       juncId=c(jens.reverse$id,jens.reverse$id))

    tabOver <- as.data.frame(findOverlaps(juncti.GR,jens.GR,ignore.strand=FALSE))
    tabOver <- cbind(tabOver,transcripts=jens[match(mcols(jens.GR[tabOver$subjectHits])[,"juncId"],jens$id),"transcript_id"])

    ## remove duplicates
    tabOver <- tabOver[!duplicated(tabOver[,c("queryHits","transcripts")]),]
    if(nrow(tabOver)>0)
    {
        setDT(tabOver)
        tabOver <- tabOver %>% group_by(queryHits) %>% summarise(transcripts = str_c(na.omit(transcripts), collapse=","))
        setDF(tabOver)
        splitReadTable[match(mcols(juncti.GR[tabOver$queryHits])[,"juncId"],splitReadTable$junID),"acceptor"] <- tabOver$transcripts
    }

    rm(tabOver,jens.GR,juncti.GR,jens.reverse)

    ## get original coordinates
    splitReadTable$start <- splitReadTable$start +1
    splitReadTable$stop <- splitReadTable$stop -1

    return(splitReadTable)
}






