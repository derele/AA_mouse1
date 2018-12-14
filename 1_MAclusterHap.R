library(MultiAmplicon)
library(pheatmap)
library(data.table)
library(DECIPHER)
library(ape)
library(pegas)

FILTER <- TRUE
newMA <- TRUE
newDeDa <- TRUE

files <- list.files(path="/SAN/Victors_playground/mouseGT/raw_fastq",
                    pattern=".fastq.gz", full.names=TRUE)

## first pool file names
Ffq.file <- files[grepl("R1_001", files)]
Rfq.file <- files[grepl("R2_001", files)]

## first pool sample names
samples <- gsub("/SAN/Victors_playground/mouseGT/raw_fastq/(.*?)_S\\d{,3}_L001_R1_001.fastq\\.gz", "\\1", Ffq.file)

filt_path <- "/SAN/Victors_playground/mouseGT/DaDaFilt"
if(!file_test("-d", filt_path)) dir.create(filt_path)

filtFs <- file.path(filt_path, paste0(samples, "_F_filt.fastq.gz"))
names(filtFs) <- samples
filtRs <- file.path(filt_path, paste0(samples, "_R_filt.fastq.gz"))
names(filtRs) <- samples

if (FILTER){
   filter.track <- lapply(seq_along(Ffq.file),  function (i) {
       filterAndTrim(Ffq.file[i], filtFs[i], Rfq.file[i], filtRs[i],
                     truncLen=c(270,270), minLen=c(270,270), 
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

names(primerF) <- as.character(ptable[,1])
names(primerR) <- as.character(ptable[,3])

files <- PairedReadFileSet(filtFs, filtRs)

primers <- PrimerPairsSet(primerF, primerR)

rownames(ptable) <- names(primers)

MA <- MultiAmplicon(primers, files)


stratfiles <- "/SAN/Victors_playground/mouseGT/stratified_files"

if(newMA) {
    if(dir.exists(stratfiles)){
        unlink(stratfiles,
               recursive=TRUE)
        }
    MA1 <- sortAmplicons(MA, filedir= stratfiles
                         starting.at=1, max.mismatch=4)
    saveRDS(MA1, file="/SAN/Zebra/MA1.Rds")
} else {
    if(!newMA){
        MA1 <- readRDS(file="/SAN/Zebra/MA1.Rds")
    } else {stop("Whant new sorting or not? Set newMA to TRUE")}
} 



pdf("figures/primers_MA_sorted.pdf", width=46)
plotAmpliconNumbers(MA1)
dev.off()


MA1 <- MA1[which(grepl("^Eim_Nuc", rownames(MA1))), ]

pdf("figures/primers_MA_ENuc.pdf", width=46)
clust <- plotAmpliconNumbers(MA1)
dev.off()

two.clusters.col <- cutree(clust$tree_col, k=2)
keep.sample <- names(two.clusters.col)[two.clusters.col==2]

MA1 <- MA1[, which(colnames(MA1)%in%keep.sample)]

pdf("figures/primers_MA_EfalClade.pdf", width=45, height=15, onefile=FALSE)
plotAmpliconNumbers(MA1)
dev.off()

## this only excludes one completely empty
MA1 <- MA1[which(rowSums(getRawCounts(MA1))>2000), ]


if(newDeDa){
    MA2 <- derepMulti(MA1, mc.cores=20)
    MA3 <- dadaMulti(MA2, Ferr=NULL, Rerr=NULL, selfConsist=TRUE,
                     multithread=FALSE, mc.cores=20, verbose=0, MAX_CONSIST=20)
    MA4.merged <- mergeMulti(MA3, justConcatenate=FALSE, mc.cores=20)
    prop.merged <- calcPropMerged(MA4.merged)
    prop.merged[is.na(prop.merged)] <- 0
    table(prop.merged>0.8)
    MA4 <- mergeMulti(MA3, justConcatenate=prop.merged<0.8,
                      mc.cores=20)
    MA5 <- makeSequenceTableMulti(MA4, mc.cores=20, orderBy="nsamples")
    foo <- removeBimeraDenovo(getSequenceTable(MA5), multithread=TRUE, method="pooled")
    MA6 <- MA5
    MA6@sequenceTableNoChime <- foo
    saveRDS(MA6, file="/SAN/Victors_playground/mouseGT/MA6.Rds")
} else {
    MA6 <- readRDS(file="/SAN/Victors_playground/mouseGT/MA6.Rds")
}

STnoC <- getSequenceTableNoChime(MA6)

## discard empty samples and ASVs for each amplicon
STnoC <- lapply(STnoC, function(x) x[rowSums(x)>0, colSums(x)>0])

pdf("figures/all_amps_heated.pdf")
lapply(STnoC, function (x) {
    pheatmap(x, show_colnames=FALSE, scale="none")
})
dev.off()


## get the most represented sequence for each amplicon and sample

mostAbundant <- lapply(STnoC, function (amp){
    apply(amp, 1, function (x){
        mostAbundant <- which(x==max(x))
        list(names(x)[mostAbundant])
    })
})

### create empty sequences for each "null allele"
mAdf <- as.data.frame(rbindlist(mostAbundant, fill=TRUE))
ISnull <- sapply(mAdf, function (x) sapply(x, is.null))
mAdf[ISnull] <- "NNN"

#### Bring Ns for each null allel on the length of the sequence
seqdf <- apply(mAdf, 1, function (x) {
    slen <- max(nchar(x))
    repl <- paste(rep("N", times=slen), collapse="")
    gsub("NNN", repl, x)
})

DNA <- apply(seqdf, 2, function (x) { 
    d <- DNAStringSet(as.character(x))
    names(d) <- colnames(mAdf)
    d
})

alnDNA1 <- lapply(DNA, AlignSeqs)

## now exclude too distant sequences (few seqs are from wrong targets
## I guess)
distDNA <- lapply(alnDNA1, DistanceMatrix, penalizeGapLetterMatches = FALSE)

## to see maximal distances
lapply(distDNA, max)
##  -> >20% difference is clearly nonsense

## so we can exclude sequences which have with more than the are in
## the minority as having huge differnces 
differ <- lapply(distDNA, function (x) colSums(x>0.2))
toExcl <- lapply(differ, function (x) which (x > length(x)/2)) 

DNA2 <- lapply(seq_along(DNA), function (i) {
    D <- DNA[[i]]
    slen <- max(width(DNA[[i]]))
    repl <- paste(rep("N", times=slen), collapse="")
    D[toExcl[[i]]] <- repl
    D
})

alnDNA2 <- lapply(DNA2, AlignSeqs)

## checking that huge distances are gone
distDNA2 <- lapply(alnDNA2, DistanceMatrix, penalizeGapLetterMatches = FALSE)

## to see maximal distances
lapply(distDNA2, max)

## to chech alignment length vs longest sequence
longSeq <- unlist(lapply(DNA, function (x) max(width(x))))
alnSeq1 <- unlist(lapply(alnDNA1, function (x) max(width(x))))
alnSeq2 <- unlist(lapply(alnDNA2, function (x) max(width(x))))
cbind(longSeq, alnSeq1, alnSeq2)

## to get the alingment partitions
cumsum(alnSeq2)

## paste the string sets togehter seq by seq
pasteDNAstrings <- function (x, y){
    DNAStringSet(paste0(as.character(x), as.character(y)))
}

oneAlnDNA <- Reduce(pasteDNAstrings, alnDNA2)
names(oneAlnDNA) <- colnames(mAdf)

writeXStringSet(oneAlnDNA, file="/SAN/Victors_playground/mouseGT/combinedEfalNuc.fasta")

dbin <- as.DNAbin(oneAlnDNA)

## library(sidier)
## mutation.network(dbin)
## mutation.network(dbinWO, network.method="zero")

hapDNAString <- DNAStringSet(apply(as.character(hap), 1, paste, collapse=""))

d <- DistanceMatrix(hapDNAString, penalizeGapLetterMatches = FALSE)

hap <- haplotype(dbin)
hap <- sort(hap, what = "label")

pdf("figures/hapMDS.pdf")
plotNetMDS(net, d)
dev.off()


assoTab <- stack(setNames(attr(hap, "index"), rownames(hap)))

assoTab$sample <- names(oneAlnDNA)[assoTab$values]

ind.hap <- table(hap=assoTab$ind, individuals=assoTab$sample)

net <- haploNet(hap, d)

pdf("figures/hapNet.pdf", width=62, height=62)
plot(net, size=attr(net, "freq"),
     scale.ratio=0.2, pie=ind.hap)
legend(-8, 0, colnames(ind.hap),
       col=rainbow(nrow(assoTab)), pch=19, ncol=2)
dev.off()

## now with all gaps deleted
## dbinWO <- del.colgapsonly(as.matrix(dbin), threshold = 0.01)


## hapWO <- haplotype(dbinWO)
## hapWO <- sort(hapWO, what = "label")

## netWO <- haploNet(hapWO)

## assoTabWO <- stack(setNames(attr(hapWO, "index"), rownames(hapWO)))

## assoTabWO$sample <- names(oneAlnDNA)[assoTabWO$values]

## ind.hapWO <- table(hap=assoTabWO$ind, individuals=assoTabWO$sample)

## pdf("figures/hapNetWO.pdf", width=62, height=62)
## plot(netWO, size=attr(netWO, "freq"),
##      scale.ratio=0.2, pie=ind.hapWO)
## legend(-8, 0, colnames(ind.hapWO),
##        col=rainbow(nrow(assoTabWO)), pch=19, ncol=2)
## dev.off()
