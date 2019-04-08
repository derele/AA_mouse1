# Alice - adapted from Emanuel's code to extract nuclear and api markers
library(MultiAmplicon)
library(pheatmap)
library(data.table)
library(DECIPHER)
library(ape)
library(pegas)
library(sidier)
library(adegenet)
library(wordcloud)
library(systemPipeR)

FILTER <- FALSE 
newMA <- FALSE
newDeDa <- FALSE

files <- list.files(path="/SAN/Victors_playground/mouseGT/raw_fastq",
                    pattern=".fastq.gz", full.names=TRUE)

## first pool file names
Ffq.file <- files[grepl("R1_001", files)]
Rfq.file <- files[grepl("R2_001", files)]

## first pool sample names
samples <- gsub("/SAN/Victors_playground/mouseGT/raw_fastq/(.*?)_S\\d{,3}_L001_R1_001.fastq\\.gz", "\\1", Ffq.file)
# -> 192 samples in total

filt_path <- "/SAN/Victors_playground/mouseGT/DaDaFilt"
if(!file_test("-d", filt_path)) dir.create(filt_path)

filtFs <- file.path(filt_path, paste0(samples, "_F_filt.fastq.gz"))
names(filtFs) <- samples
filtRs <- file.path(filt_path, paste0(samples, "_R_filt.fastq.gz"))
names(filtRs) <- samples

if (FILTER){
  filter.track <- lapply(seq_along(Ffq.file),  function (i) {
    filterAndTrim(Ffq.file[i], filtFs[i], Rfq.file[i], filtRs[i],
                  truncLen=c(285,285), minLen=c(285,285), 
                  maxN=0, maxEE=2, truncQ=2, 
                  compress=TRUE, verbose=TRUE)
  })
  saveRDS(filter.track, file="/SAN/Victors_playground/mouseGT/filter.Rds")
} else {
  filter.track <- readRDS(file="/SAN/Victors_playground/mouseGT/filter.Rds")
}

filter <- do.call(rbind, filter.track)
colSums(filter)[2]/colSums(filter)[1]
## low proportion of data overall cept (but okay for this purpose!)

file.ptable <- "/SAN/Victors_playground/mouseGT/primer.file.csv"

ptable <- read.csv(file.ptable, sep=",", header=TRUE)

primerF <- gsub(" ", "", ptable[,2])
primerR <- gsub(" ", "", ptable[,4])

# correct typo space!
names(primerF) <- gsub(" ", "", as.character(ptable[,1]))
names(primerR) <- gsub(" ", "", as.character(ptable[,3]))

files <- PairedReadFileSet(filtFs, filtRs)

primers <- PrimerPairsSet(primerF, primerR)

rownames(ptable) <- names(primers)

MA <- MultiAmplicon(primers, files)

stratfiles <- "/SAN/Alices_sandpit/AA_EimeriaMarkers/data/stratified_files/"

if(newMA) {
  if(dir.exists(stratfiles)){
    unlink(stratfiles,
           recursive=TRUE)
  }
  MA1 <- sortAmplicons(MA, filedir=stratfiles,
                       starting.at=1, max.mismatch=4)
  saveRDS(MA1, file="/SAN/Alices_sandpit/AA_EimeriaMarkers/data/MA1.Rds")
  pdf("/SAN/Alices_sandpit/AA_EimeriaMarkers/data/figures/primers_MA_sorted.pdf", width=46)
  plotAmpliconNumbers(MA1)
  dev.off()
} else {
  if(!newMA){
    MA1 <- readRDS(file="/SAN/Alices_sandpit/AA_EimeriaMarkers/data/MA1.Rds")
  } else {stop("Whant new sorting or not? Set newMA to TRUE")}
} 

## MA1 contains the full filtrated data. From then on, sep api / nucl
## Apicoplast
MA1.api <- MA1[which(grepl("^Eim_Ap", rownames(MA1))), ]
## Nucl
MA1.nucl <- MA1[which(grepl("^Eim_Nuc", rownames(MA1))), ]

## Nucl
clust.nucl <- plotAmpliconNumbers(MA1.nucl)

# cut trees in 2 clusters 
## Nucl
two.clusters.col.nucl <- cutree(clust.nucl$tree_col, k=2)
keep.sample.nucl <- names(two.clusters.col.nucl)[two.clusters.col.nucl==2]
MA1.nucl <- MA1.nucl[, which(colnames(MA1.nucl)%in%keep.sample.nucl)]

plotAmpliconNumbers(MA1.nucl)

# newDeDa was set to false
if(newDeDa){
  MA2.nucl <- derepMulti(MA1.nucl, mc.cores=20)
  MA3.nucl <- dadaMulti(MA2.nucl, Ferr=NULL, Rerr=NULL, selfConsist=TRUE,
                   multithread=FALSE, mc.cores=20, verbose=0, MAX_CONSIST=20)
  MA4.nucl <- mergeMulti(MA3.nucl, mc.cores=20)
  MA5.nucl <- makeSequenceTableMulti(MA4.nucl, mc.cores=20, orderBy="nsamples")
  MA6.nucl <- removeChimeraMulti(MA5.nucl)
  saveRDS(MA6.nucl, file="/SAN/Alices_sandpit/AA_EimeriaMarkers/data/MA6.nucl.Rds")
} else {
  MA6.nucl <- readRDS(file="/SAN/Alices_sandpit/AA_EimeriaMarkers/data/MA6.nucl.Rds")
}

## this only excludes one completely empty
## subsetting is broken in package !!! Jan 5th.. 
## MA6 <- MA6[which(rowSums(getRawCounts(MA6))>2000), ]

STnoC.nucl <- getSequenceTableNoChime(MA6.nucl)

## discard empty samples and ASVs for each amplicon
STnoC.nucl <- lapply(STnoC.nucl, function(x) x[rowSums(x)>0, colSums(x)>0])

STnoC.nucl <- STnoC[unlist(lapply(STnoC.nucl, function (x) all(dim(x) >1)))]
# 
# pdf("figures/all_amps_heated.pdf")
# lapply(STnoC, function (x) {
#   pheatmap(x, show_colnames=FALSE, scale="none")
# })
# dev.off()

## get the most represented sequence for each amplicon and sample
DNA.nucl <- lapply(STnoC.nucl, function (amp){
  mostAbSeq <- apply(amp, 1, function (x){
    mostAbundant <- which(x==max(x))
    names(x)[mostAbundant]
  })
  DNAStringSet(mostAbSeq)
})

## NB remove empty primer
lengthMarkers.nucl <- lapply(DNA.nucl, function (x) length(DNAStringSet(gsub("-", "", x))))
DNA.nucl <- DNA.nucl[lengthMarkers.nucl != 0]

## now we can remove the Ns and just let the alignment put gaps
alnDNA1.nucl <- lapply(DNA.nucl, function (x) AlignSeqs(x))

## now exclude too distant sequences (few seqs are from wrong targets
## I guess)
DNA2.nucl <- lapply(alnDNA1.nucl, function(x){
  distDNA <- DistanceMatrix(x, penalizeGapLetterMatches = FALSE)
  d <- colSums(distDNA>0.2, na.rm=TRUE)
  ## so we can exclude sequences which have with more than the are in
  ## the minority as having huge differnces 
  toKeep <- which (d < (length(x)/2)) 
  x[toKeep]
})

## realign removeing the gaps before
alnDNA2.nucl <- lapply(DNA2.nucl, function (x) AlignSeqs(DNAStringSet(gsub("-", "", x))))

## checking that huge distances are gone
distDNA2.nucl <- lapply(alnDNA2.nucl, DistanceMatrix, penalizeGapLetterMatches = FALSE)

## to see maximal distances
summary(unlist(lapply(distDNA2.nucl, max, na.rm=TRUE)))

## to chech alignment length vs longest sequence
longSeq.nucl <- unlist(lapply(DNA.nucl, function (x) max(width(x))))
alnSeq1.nucl <- unlist(lapply(alnDNA1.nucl, function (x) max(width(x))))
alnSeq2.nucl <- unlist(lapply(alnDNA2.nucl, function (x) max(width(x))))
cbind(longSeq.nucl, alnSeq1.nucl, alnSeq2.nucl)

## remove the worst markers producing hardly any sequence for lots of
## samples

## remove duplicate technical replicate sequences for all samples
## return NNN sequences to have all alignments have the same sequences

all.reps.nucl <- unique(unlist(lapply(alnDNA2.nucl, names)))
all.samples.nucl <- as.factor(gsub("-chip\\d", "", all.reps.nucl))

Aln.nucl <- lapply(alnDNA2.nucl, function (aln){
  alnWidth <- unique(width(aln))
  sampleNames <- gsub("-chip\\d", "", names(aln))
  ## could probably also work with an empty DNAstring("")
  badString <- DNAString(paste(rep("N", times=alnWidth), collapse=""))
  marker <- sapply(levels(all.samples.nucl), function (x){
    pair <- which(sampleNames%in%x)
    DNA <- aln[pair]
    ## ## if it's only one sequence return it
    if(length(DNA)<1){
      return(badString)
    }
    if(all(!grepl("[A|G|C|T]", DNA))){return(badString)}
    if(length(DNA)<2){return(DNA[[1]])}
    ## ## if both are gap only return one of them (the first)        
    ## PROBLEMATIC: What if the two technical reps don't agree?
    ## repD <- DistanceMatrix(DNA)
    ## if(any(repD>0)){
    ##     warning(paste("replicate dist >0:", repD[repD>0][1]))
    ##     return(DNAString(paste(rep("N", times=alnWidth),
    ##                            collapse="")))
    ## }
    else {
      return(DNA[[1]])
    }
  })
  DNAStringSet(marker) 
})

## throw out markers with less then optimal number of sequences
table(unlist(lapply(Aln.nucl, function (x) length(grep("[A|G|C|T]", x)))))

Aln.nucl <- Aln.nucl[unlist(lapply(Aln.nucl, function (x) length(grep("[A|G|C|T]", x))))>13]

## to get the alingment partitions
partitionMarkers.nucl <- cumsum(unlist(lapply(Aln.nucl, function (x) max(width(x)))))

## paste the string sets togehter seq by seq
pasteDNAstrings <- function (x, y){
  DNAStringSet(paste0(as.character(x), as.character(y)))
}

oneAln.nucl <- Reduce(pasteDNAstrings, Aln.nucl)
names(oneAln.nucl) <- names(Aln.nucl[[1]])

# And save
pm.nucl <- c(start = 0, partitionMarkers.nucl)
for (marker in 1:(length(pm.nucl)-1)){
  x <- list()
  for (sample in 1:length(oneAln.nucl)){
    seq <- oneAln.nucl[[sample]][pm.nucl[marker] : pm.nucl[marker+1]]
    x <- c(x, paste0(">", names(oneAln.nucl[sample]), "\n", seq, "\n"))
  }
  write.table(unlist(x),
              paste0("/SAN/Alices_sandpit/AA_EimeriaMarkers/data/seqNucl/seq_", names(pm.nucl)[marker+1], ".fasta"),
              quote = F, row.names = F, eol = "", col.names = F)
}

## 33 samples now???

#################################
## And then, on apicoplast markers, same 32 samples, extract and align
samplesGoodForNucl <- names(oneAln.nucl)

MA1.api <- MA1.api[,grep(paste(samplesGoodForNucl, collapse="|"), colnames(MA1.api))]

plotAmpliconNumbers(MA1.api)

newDeDa2 = FALSE
if(newDeDa2){
  MA2.api <- derepMulti(MA1.api, mc.cores=20)
  MA3.api <- dadaMulti(MA2.api, Ferr=NULL, Rerr=NULL, selfConsist=TRUE,
                   multithread=FALSE, mc.cores=20, verbose=0, MAX_CONSIST=20)
  MA4.api <- mergeMulti(MA3.api, mc.cores=20)
  MA5.api <- makeSequenceTableMulti(MA4.api, mc.cores=20, orderBy="nsamples")
  MA6.api <- removeChimeraMulti(MA5.api)
  saveRDS(MA6.api, file="/SAN/Alices_sandpit/AA_EimeriaMarkers/data/MA6.api.Rds")
} else {
  MA6.api <- readRDS(file="/SAN/Alices_sandpit/AA_EimeriaMarkers/data/MA6.api.Rds")
}

STnoC.api <- getSequenceTableNoChime(MA6.api)

## discard empty samples and ASVs for each amplicon
STnoC.api <- lapply(STnoC.api, function(x) x[rowSums(x)>0, colSums(x)>0])

STnoC.api <- STnoC.api[unlist(lapply(STnoC.api, function (x) all(dim(x) >1)))]

## get the most represented sequence for each amplicon and sample
DNA.api <- lapply(STnoC.api, function (amp){
  mostAbSeq <- apply(amp, 1, function (x){
    mostAbundant <- which(x==max(x))
    names(x)[mostAbundant]
  })
  DNAStringSet(mostAbSeq)
})

## now we can remove the Ns and just let the alignment put gaps
alnDNA1.api <- lapply(DNA.api, function (x) AlignSeqs(x))

## now exclude too distant sequences (few seqs are from wrong targets
## I guess)

DNA2.api <- lapply(alnDNA1.api, function(x){
  distDNA <- DistanceMatrix(x, penalizeGapLetterMatches = FALSE)
  d <- colSums(distDNA>0.2, na.rm=TRUE)
  ## so we can exclude sequences which have with more than the are in
  ## the minority as having huge differnces 
  toKeep <- which (d < (length(x)/2)) 
  x[toKeep]
})

## NB remove empty primer
lengthMarkers.api <- lapply(DNA2.api, function (x) length(DNAStringSet(gsub("-", "", x))))
DNA2.api <- DNA2.api[lengthMarkers.api != 0]

## realign removeing the gaps before
alnDNA2.api <- lapply(DNA2.api, function (x) AlignSeqs(DNAStringSet(gsub("-", "", x))))

## checking that huge distances are gone
distDNA2.api <- lapply(alnDNA2.api, DistanceMatrix, penalizeGapLetterMatches = FALSE)

## to see maximal distances
summary(unlist(lapply(distDNA2.api, max, na.rm=TRUE)))

## to chech alignment length vs longest sequence
longSeq.api <- unlist(lapply(DNA.api, function (x) max(width(x))))
alnSeq1.api <- unlist(lapply(alnDNA1.api, function (x) max(width(x))))
alnSeq2.api <- unlist(lapply(alnDNA2.api, function (x) max(width(x))))
cbind(longSeq.api, alnSeq1.api, alnSeq2.api)

## remove the worst markers producing hardly any sequence for lots of
## samples

## remove duplicate technical replicate sequences for all samples
## return NNN sequences to have all alignments have the same sequences

all.reps.api <- unique(unlist(lapply(alnDNA2.api, names)))
all.samples.api <- as.factor(gsub("-chip\\d", "", all.reps.api))

Aln.api <- lapply(alnDNA2.api, function (aln){
  alnWidth <- unique(width(aln))
  sampleNames <- gsub("-chip\\d", "", names(aln))
  ## could probably also work with an empty DNAstring("")
  badString <- DNAString(paste(rep("N", times=alnWidth), collapse=""))
  marker <- sapply(levels(all.samples.api), function (x){
    pair <- which(sampleNames%in%x)
    DNA <- aln[pair]
    ## ## if it's only one sequence return it
    if(length(DNA)<1){
      return(badString)
    }
    if(all(!grepl("[A|G|C|T]", DNA))){return(badString)}
    if(length(DNA)<2){return(DNA[[1]])}
    ## ## if both are gap only return one of them (the first)        
    ## PROBLEMATIC: What if the two technical reps don't agree?
    ## repD <- DistanceMatrix(DNA)
    ## if(any(repD>0)){
    ##     warning(paste("replicate dist >0:", repD[repD>0][1]))
    ##     return(DNAString(paste(rep("N", times=alnWidth),
    ##                            collapse="")))
    ## }
    else {
      return(DNA[[1]])
    }
  })
  DNAStringSet(marker) 
})

## throw out markers with less then optimal number of sequences
table(unlist(lapply(Aln.api, function (x) length(grep("[A|G|C|T]", x)))))

Aln.api <- Aln.api[unlist(lapply(Aln.api, function (x) length(grep("[A|G|C|T]", x))))>13]

## to get the alingment partitions
partitionMarkers.api <- cumsum(unlist(lapply(Aln.api, function (x) max(width(x)))))

## paste the string sets togehter seq by seq
pasteDNAstrings <- function (x, y){
  DNAStringSet(paste0(as.character(x), as.character(y)))
}

oneAln.api <- Reduce(pasteDNAstrings, Aln.api)
names(oneAln.api) <- names(Aln.api[[1]])

pm.api <- c(start = 0, partitionMarkers.api)

## For all markers:
for (marker in 1:(length(pm.api)-1)){
  x <- list()
  for (sample in 1:length(oneAln.api)){
    seq <- oneAln.api[[sample]][pm.api[marker] : pm.api[marker+1]]
    x <- c(x, paste0(">", names(oneAln.api[sample]), "\n", seq, "\n"))
  }
  write.table(unlist(x),
              paste0("/SAN/Alices_sandpit/AA_EimeriaMarkers/data/seqApi/seq_", names(pm.api)[marker+1], ".fasta"),
              quote = F, row.names = F, eol = "", col.names = F)
}

# to complete with Victor
write.csv(data.frame(marker = names(partitionMarkers.api), pos = partitionMarkers.api, frame = NA), file = "/SAN/Alices_sandpit/AA_EimeriaMarkers/combinedEfalApiFramePartitions.csv", row.names = F, quote = F)
