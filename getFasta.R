## Alice April 2019
## after extract markers (R) and blast (bash)
library(ShortRead)

markers <- gsub(".fasta.blastn", "",
                gsub("/SAN/Alices_sandpit/AA_mouse1/git/AA_mouse1/blast/Efal/", "",
                     list.files("/SAN/Alices_sandpit/AA_mouse1/git/AA_mouse1/blast/Efal/")))

markers.api <- markers[grep("Eim_Ap", markers)]
markers.nucl <- markers[grep("Eim_Nuc", markers)]

# Get the raw sequences:
dfFull.nucl <- data.frame()
for (marker in markers.nucl){
  df <- readFasta(paste0("/SAN/Alices_sandpit/AA_mouse1/data/seqNucl/",
                         marker, ".fasta"))
  df <- data.frame(FullSampleSeq = df@sread, marker = marker, sample = df@id)
  dfFull.nucl <- rbind(dfFull.nucl, df)
}

dfFull.api <- data.frame()
for (marker in markers.api){
  df <- readFasta(paste0("/SAN/Alices_sandpit/AA_mouse1/data/seqApi/",
                         marker, ".fasta"))
  df <- data.frame(FullSampleSeq = df@sread, marker = marker, sample = df@id)
  dfFull.api <- rbind(dfFull.api, df)
}

# Get the blast results:
getBlastResults <- function(isolate){
  blastFull <- data.frame()
  for (marker in markers){
    f <- paste0("/SAN/Alices_sandpit/AA_mouse1/git/AA_mouse1/blast/", isolate, "/",
                marker, ".fasta.blastn")
    if (file.size(f) > 0){
      blast <- read.csv(f, header = F)
      blast$marker <- marker
      blastFull <- rbind(blastFull, blast)
    }
  }
  nameMyBlast <- function(blastdf){
    names(blastdf) <-  c("sample","saccver", "pident", "length", "mismatch", "gapopen", "qstart", "qend",
                         "sstart", "send", "evalue", "bitscore", "staxid", isolate, "queryseq", "qframe", "sframe", "marker")
    return(blastdf)
  }
  blastFull <- nameMyBlast(blastFull)
  blastFull.Api <- blastFull[grep("Eim_Ap", blastFull$marker),]
  blastFull.Nuc <- blastFull[grep("Eim_Nuc", blastFull$marker),]
  return(list(blastFull = blastFull, blastFull.Api = blastFull.Api, blastFull.Nuc = blastFull.Nuc))
}

blastList_E139 <- getBlastResults(isolate = "E139")
blastList_E88 <- getBlastResults(isolate = "E88")
# errors for these 2, crapy genomes :( I need to go back there
blastList_E64 <- getBlastResults(isolate = "E64")
blastList_Efal <- getBlastResults(isolate = "Efal")
blastList_Ever <- getBlastResults(isolate = "Ever")

# Merge both

# Efal
all.Nuc.Efal <- merge(blastList_Efal$blastFull.Nuc, dfFull.nucl)
length(unique(all.Nuc.Efal$marker)) # ok
all.Api.Efal <- merge(blastList_Efal$blastFull.Api, dfFull.api)
length(unique(all.Api.Efal$marker)) # ok

# Ever
all.Nuc.Ever <- merge(blastList_Ever$blastFull.Nuc, dfFull.nucl)
length(unique(all.Nuc.Ever$marker)) # ok
all.Api.Ever <- merge(blastList_Ever$blastFull.Api, dfFull.api)
length(unique(all.Api.Ever$marker)) # ok

# E64
all.Nuc.E64 <- merge(blastList_E64$blastFull.Nuc, dfFull.nucl)
length(unique(all.Nuc.E64$marker)) # only 20 were found by blast
all.Api.E64 <- merge(blastList_E64$blastFull.Api, dfFull.api)
length(unique(all.Api.E64$marker)) # ok

# select the reference sequences to keep
## to put them in the order we need partition files
# (also found for nucl in /SAN/Victors_playground/mouseGT/combinedEfalNucFramePartitions.csv)
# incomplete for api, was done during extractMarkers.R by Alice
fororderNuc <- read.csv("/SAN/Alices_sandpit/AA_mouse1/git/AA_mouse1/combinedEfalNucFramePartitions.csv")
fororderApi <- read.csv("/SAN/Alices_sandpit/AA_mouse1/git/AA_mouse1/combinedEfalApiFramePartitions.csv")

names(fororderNuc)[1] <- "marker"
fororderNuc$rank <- 1:nrow(fororderNuc)
fororderNuc$marker <- as.character(fororderNuc$marker)
names(fororderApi)[1] <- "marker"
fororderApi$rank <- 1:nrow(fororderApi)
fororderApi$marker <- as.character(fororderApi$marker)

library(dplyr)

getRefSeq2Keep <- function(isolate, genome){
  fororder = get(paste0("fororder", genome))
  all = get(paste0("all.", genome, ".", isolate))
  final = data.frame(all %>%
                       distinct(marker, .keep_all = TRUE) %>%
                       group_by(marker) %>% # per marker
                       dplyr::filter(evalue == min(evalue)) %>% 
                       dplyr::filter(pident == max(pident)) %>% 
                       dplyr::filter(bitscore == max(bitscore)))
  final$marker = gsub("seq_", "", final$marker)
  final <- merge(final, fororder) # to order
  final <- final[order(final$rank), ]
  
  longSeq = paste0(">", genome, isolate, "allMarkers\n",
             as.character(apply(final[isolate], 2, function(row) paste(row, collapse = ""))))
  write.table(x = longSeq,
              paste0("/SAN/Alices_sandpit/AA_mouse1/git/AA_mouse1/markersByRef/", "markersFrom", isolate, ".", genome, ".fasta"),
              row.names =  F, quote = F, col.names = F)
  
  return(list(final= final, longSeq = longSeq))
}

final.nucl.Efal <- getRefSeq2Keep("Efal", "Nuc")
final.api.Efal <- getRefSeq2Keep("Efal", "Api")
final.nucl.Ever <- getRefSeq2Keep("Ever", "Nuc")
final.api.Ever <- getRefSeq2Keep("Ever", "Api")
final.nucl.E64 <- getRefSeq2Keep("E64", "Nuc")
final.api.E64 <- getRefSeq2Keep("E64", "Api")
