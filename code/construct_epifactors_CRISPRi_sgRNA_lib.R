#Project: iPS CRISPRi
#construct_epifactors_CRISPRi_sgRNA_lib.R
#Objective: generate a CRISPRi sgRNA library targeting all known epigenetic factors in the human genome
#Mohamad Najia

library(ggplot2)
library(stringr)
library(dplyr)
library(data.table)


### Initialization ###

#import library annotation file 
epifactors_main <- fread(file = "./epifactors_sgRNA_library/EpiGenes_main_v1.7.3.tsv", data.table = FALSE)
rownames(epifactors_main) <- epifactors_main$HGNC_symbol

#import genome-wide CRISPRi libraries 
dolcetto <- fread(file = "./epifactors_sgRNA_library/broadgpp-dolcetto-targets-seta.txt", data.table = FALSE)
hCRISPRi <- fread(file = "./epifactors_sgRNA_library/hCRISPRiv2.1_libraries.tsv", data.table = FALSE)

#specify 5' and 3' constant sequences for the oligo pool 
constant5 <- "AAGCAGTGGTATCAACGCAGAGTcgtctcgcacc"
constant3 <- "gtttagagacgCCTTTGGGTAAGCACACGTC"


### Construct CRISPRi sgRNA Library ###

#first, extract sgRNAs from the Dolcetto library where the gene matches the GeneID in the epifactors annotation file 
dmatches <- dolcetto[dolcetto$`Annotated Gene ID` %in% epifactors_main$GeneID,]
colnames(dmatches) <- c("sgRNA_Sequence", "Gene_Symbol", "GeneID")
dmatches$sgRNA_Type <- "gene_target"
dmatches$Source <- "Dolcetto"
dmatches$Library_Subset <- "A"
dmatches$Transcript <- ""

#next, extract sgRNAs from the hCRISPRi_v2.1 library for the genes that were not found in the Dolcetto library
unmatched <- epifactors_main[epifactors_main$GeneID %in% setdiff(epifactors_main$GeneID, dmatches$GeneID),2]
hmatches <- hCRISPRi[hCRISPRi$gene %in% unmatched,]
hmatches <- hmatches[hmatches$`Sublibrary half` == "Top5",c("gene","transcript","protospacer sequence","Sublibrary half")]
colnames(hmatches) <- c("Gene_Symbol", "Transcript", "sgRNA_Sequence", "Library_Subset")
hmatches$sgRNA_Type <- "gene_target"
hmatches$Source <- "hCRISPRi-v2.1"
hmatches$GeneID <- epifactors_main[hmatches$Gene_Symbol, "GeneID"]

#include non-targeting control sgRNAs 
nontarget_sgRNA <- fread(file = "./epifactors_sgRNA_library/non_targeting_sgRNAs_dolcetto.tsv", header = FALSE, data.table = FALSE)
colnames(nontarget_sgRNA) <- c("Gene_Symbol", "sgRNA_Sequence")
nontarget_sgRNA$GeneID <- nontarget_sgRNA$Gene_Symbol
nontarget_sgRNA$sgRNA_Type <- "non_targeting_control"
nontarget_sgRNA$Source <- "Dolcetto"
nontarget_sgRNA$Library_Subset <- "A"
nontarget_sgRNA$Transcript <- ""

#merge sgRNAs into a single library 
epifactors_lib <- rbind(dmatches, hmatches, nontarget_sgRNA)

#identify sgRNA sequences that contain BsmBI sites and remove from the final library
bsmbi <- "CGTCTC"
bsmbirc <- "GAGACG"
inds <- grepl(bsmbi, epifactors_lib$sgRNA_Sequence, ignore.case = TRUE) | grepl(bsmbirc, epifactors_lib$sgRNA_Sequence, ignore.case = TRUE)
epifactors_lib <- epifactors_lib[!inds,]
epifactors_lib$Contains_BsmBI <- "No"

#create sgRNA sequences for the oligo pool
epifactors_lib$Synthesized_sgRNA <- epifactors_lib$sgRNA_Sequence
notG <- substr(epifactors_lib$sgRNA_Sequence, 1,1) != "G"
epifactors_lib$Synthesized_sgRNA[notG] <- paste("G", epifactors_lib$Synthesized_sgRNA[notG], sep = "")
epifactors_lib$Oligo_Sequence <- paste(constant5, epifactors_lib$Synthesized_sgRNA, constant3, sep = "")
epifactors_lib$sgRNA_Length <- nchar(epifactors_lib$Synthesized_sgRNA)

#tidying up 
epifactors_lib <- epifactors_lib[order(epifactors_lib$Gene_Symbol),]
rownames(epifactors_lib) <- 1:dim(epifactors_lib)[1]
epifactors_lib$sgRNA_ID <- paste("EpiKDlib", epifactors_lib$Gene_Symbol, rownames(epifactors_lib), sep = "_")
epifactors_lib <- epifactors_lib[, c("Gene_Symbol", "GeneID", "Transcript", "sgRNA_ID", "sgRNA_Sequence", "Synthesized_sgRNA", "Oligo_Sequence", "sgRNA_Length", "sgRNA_Type", "Source", "Library_Subset", "Contains_BsmBI")]

#write the sgRNA library to disk
write.table(epifactors_lib, file = "./epifactors_sgRNA_library/EpiFactors_CRISPRi_lib_v3.tsv", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)




