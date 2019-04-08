## Alice April 2019
## after extract markers (R) and blast (bash)
library(ShortRead)

markers <- gsub(".fasta.blastn", "",
                gsub("/SAN/Alices_sandpit/AA_EimeriaMarkers/blast/Efal/", "",
                     list.files("/SAN/Alices_sandpit/AA_EimeriaMarkers/blast/Efal/")))

markers.api <- markers[grep("Eim_Ap", markers)]
markers.nucl <- markers[grep("Eim_Nuc", markers)]

# Get the blast results:
blastFullEfal <- data.frame()
for (marker in markers){
  blast <- read.csv(paste0("/SAN/Alices_sandpit/AA_EimeriaMarkers/blast/Efal/",
                           marker, ".fasta.blastn"),
                    header = F)
  blast$marker <- marker
  blastFullEfal <- rbind(blastFullEfal, blast)
}

blastFullEfal.Api <- blastFullEfal[grep("Eim_Ap", blastFullEfal$marker),]
blastFullEfal.Nuc <- blastFullEfal[grep("Eim_Nuc", blastFullEfal$marker),]
length(unique(blastFullEfal$marker)) # ok

blastFullEver <- data.frame()
for (marker in markers){
    blast <- read.csv(paste0("/SAN/Alices_sandpit/AA_EimeriaMarkers/blast/Ever/",
                           marker, ".fasta.blastn"),
                    header = F)
  blast$marker <- marker
  blastFullEver <- rbind(blastFullEver, blast)
}

blastFullEver.Api <- blastFullEver[grep("Eim_Ap", blastFullEver$marker),]
blastFullEver.Nuc <- blastFullEver[grep("Eim_Nuc", blastFullEver$marker),]
length(unique(blastFullEver$marker)) # ok

nameMyBlast <- function(blastdf, ref){
  names(blastdf) <-  c("sample","saccver", "pident", "length", "mismatch", "gapopen", "qstart", "qend",
                       "sstart", "send", "evalue", "bitscore", "staxid", ref, "queryseq", "qframe", "sframe", "marker")
  return(blastdf)
}

blastFullEfal.Api <- nameMyBlast(blastFullEfal.Api, "Efalseq")
blastFullEfal.Nuc <- nameMyBlast(blastFullEfal.Nuc, "Efalseq")
blastFullEver.Api <- nameMyBlast(blastFullEver.Api, "Everseq")
blastFullEver.Nuc <- nameMyBlast(blastFullEver.Nuc, "Everseq")

# Get the raw sequences:
dfFull.nucl <- data.frame()
for (marker in markers.nucl){
  df <- readFasta(paste0("/SAN/Alices_sandpit/AA_EimeriaMarkers/data/seqNucl/",
                           marker, ".fasta"))
  df <- data.frame(FullSampleSeq = df@sread, marker = marker, sample = df@id)
  dfFull.nucl <- rbind(dfFull.nucl, df)
}

dfFull.api <- data.frame()
for (marker in markers.api){
  df <- readFasta(paste0("/SAN/Alices_sandpit/AA_EimeriaMarkers/data/seqApi/",
                         marker, ".fasta"))
  df <- data.frame(FullSampleSeq = df@sread, marker = marker, sample = df@id)
  dfFull.api <- rbind(dfFull.api, df)
}

# merge both
all.nuc.Efal <- merge(blastFullEfal.Nuc, dfFull.nucl)
length(unique(all.nuc.Efal$marker)) # ok
all.api.Efal <- merge(blastFullEfal.Api, dfFull.api)
length(unique(all.api.Efal$marker)) # ok
all.nuc.Ever <- merge(blastFullEver.Nuc, dfFull.nucl)
length(unique(all.nuc.Ever$marker)) # ok
all.api.Ever <- merge(blastFullEver.Api, dfFull.api)
length(unique(all.api.Ever$marker)) # ok

# select the reference sequences to keep
library(dplyr)

final.nucl.Efal <- all.nuc.Efal %>%
  distinct(Efalseq, marker, .keep_all = TRUE) %>%
  group_by(marker) %>% # per marker
  dplyr::filter(evalue == min(evalue)) %>% 
  dplyr::filter(pident == max(pident)) %>% 
  dplyr::filter(bitscore == max(bitscore)) 
  
final.api.Efal <- all.api.Efal %>%
  distinct(Efalseq, marker, .keep_all = TRUE) %>%
  group_by(marker) %>% # per marker
  dplyr::filter(evalue == min(evalue)) %>% 
  dplyr::filter(pident == max(pident)) %>% 
  dplyr::filter(bitscore == max(bitscore)) 

final.nucl.Ever <- all.nuc.Ever %>%
  distinct(Everseq, marker, .keep_all = TRUE) %>%
  group_by(marker) %>% # per marker
  dplyr::filter(evalue == min(evalue)) %>% 
  dplyr::filter(pident == max(pident)) %>% 
  dplyr::filter(bitscore == max(bitscore)) 

final.api.Ever <- all.api.Ever %>%
  distinct(Everseq, marker, .keep_all = TRUE) %>%
  group_by(marker) %>% # per marker
  dplyr::filter(evalue == min(evalue)) %>% 
  dplyr::filter(pident == max(pident)) %>% 
  dplyr::filter(bitscore == max(bitscore)) 

## to put them in the order we need partition files
# (also found for nucl in /SAN/Victors_playground/mouseGT/combinedEfalNucFramePartitions.csv)
# incomplete for api, was done during extractMarkers.R by Alice

fororderNucl <- read.csv("/SAN/Alices_sandpit/AA_EimeriaMarkers/combinedEfalNucFramePartitions.csv")
fororderApi <- read.csv("/SAN/Alices_sandpit/AA_EimeriaMarkers/combinedEfalApiFramePartitions.csv")

names(fororderNucl)[1] <- "marker"
fororderNucl$rank <- 1:nrow(fororderNucl)
fororderNucl$marker <- as.character(fororderNucl$marker)
names(fororderApi)[1] <- "marker"
fororderApi$rank <- 1:nrow(fororderApi)
fororderApi$marker <- as.character(fororderApi$marker)

final.nucl.Efal$marker <- gsub("seq_", "", final.nucl.Efal$marker)
final.nucl.Efal <- merge(final.nucl.Efal, fororderNucl)
final.nucl.Efal <- final.nucl.Efal[order(final.nucl.Efal$rank), ]
X <- paste0(">allMarkersNuclFromEfal\n", 
            as.character(apply(final.nucl.Efal["Efalseq"], 2, function(row) paste(row, collapse = ""))))
write.table(x = X, 
            "/SAN/Alices_sandpit/AA_EimeriaMarkers/markersFromEfalNucl.fasta", 
            row.names =  F, quote = F, col.names = F)

final.api.Efal$marker <- gsub("seq_", "", final.api.Efal$marker)
final.api.Efal <- merge(final.api.Efal, fororderApi)
final.api.Efal <- final.api.Efal[order(final.api.Efal$rank), ]
X <- paste0(">allMarkersApiFromEfal\n", 
            as.character(apply(final.api.Efal["Efalseq"], 2, function(row) paste(row, collapse = ""))))
write.table(x = X, 
            "/SAN/Alices_sandpit/AA_EimeriaMarkers/markersFromEfalApi.fasta", 
            row.names =  F, quote = F, col.names = F)

final.nucl.Ever$marker <- gsub("seq_", "", final.nucl.Ever$marker)
final.nucl.Ever <- merge(final.nucl.Ever, fororderNucl)
final.nucl.Ever <- final.nucl.Ever[order(final.nucl.Ever$rank), ]
X <- paste0(">allMarkersNuclFromEver\n", 
            as.character(apply(final.nucl.Ever["Everseq"], 2, function(row) paste(row, collapse = ""))))
write.table(x = X, 
            "/SAN/Alices_sandpit/AA_EimeriaMarkers/markersFromEverNucl.fasta", 
            row.names =  F, quote = F, col.names = F)

final.api.Ever$marker <- gsub("seq_", "", final.api.Ever$marker)
final.api.Ever <- merge(final.api.Ever, fororderApi)
final.api.Ever <- final.api.Ever[order(final.api.Ever$rank), ]
X <- paste0(">allMarkersApiFromEver\n", 
            as.character(apply(final.api.Ever["Everseq"], 2, function(row) paste(row, collapse = ""))))
write.table(x = X, 
            "/SAN/Alices_sandpit/AA_EimeriaMarkers/markersFromEverApi.fasta", 
            row.names =  F, quote = F, col.names = F)