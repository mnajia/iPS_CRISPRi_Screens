#CRISPR Library Analysis
#Quantify sgRNA counts from NGS reads
#sgRNA_count_v4.1.R
#Mohamad Najia

library(ShortRead)
library(Biostrings)
library(stringr)


#parse input arguments 
args <- commandArgs(TRUE)
fq_file <- args[1]
sgrnas_lib_file <- args[2]
output_file <- args[3]

#import sgRNA library
sgrnas_lib <- read.table(sgrnas_lib_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
sgrnas_lib$sgRNA_Sequence <- toupper(sgrnas_lib$sgRNA_Sequence)
sgRNAs_list <- as.list(sgrnas_lib$sgRNA_Sequence)

#import fastq files 
fq <- readFastq(fq_file)
reads <- sread(fq)
reads <- as.character(reads)
total_reads <- length(reads)
print(paste("Total reads: ", total_reads, sep = ""))

#find "CACCG" TSS motif at position 20 of the read 
nums <- str_locate(reads, "CACCG")
reads <- reads[!is.na(nums[,1])]
nums <- nums[ !is.na(nums[,1]) ,]
inds <- nums[,1] == 20
print(paste0("Reads with detectable CACCG TSS motif: ", sum(inds), " (", sum(inds)/total_reads*100, "%)"))
nums <- nums[inds,]
reads <- reads[inds]

#extract 21 bps following the motif 
reads <- substr(reads, nums[,2], nums[,2]+20)

#map reads to sgRNAs based on exact matches only 
print("Mapping reads to sgRNA library with exact matches")
sample_counts <- sgrnas_lib[,c("sgRNA_ID","sgRNA_Sequence")]
truths <- rep(FALSE, length(reads))
counts_exact <- lapply(sgRNAs_list, function(x) {
  hits <- str_detect(reads, x)
  truths <<- truths | hits
  return( sum(hits) )
})

reads_map_exact <- sum(truths)
print(paste0("Reads with exact sgRNA matches: ", reads_map_exact, " (", reads_map_exact/total_reads*100, "%)"))

#write sample counts to disk 
sample_counts$sample <- unlist(counts_exact)
saveRDS(sample_counts, file = output_file)















