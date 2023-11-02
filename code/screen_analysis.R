#Project: iPS CRISPRi 
#screen_analysis.R
#Objective: Analyze iPS SSEA-5 differentiation screen of epigenetic factors
#Mohamad Najia

library(Biostrings)
library(parallel)
library(data.table)
library(LSD)
library(ggplot2)
library(ggExtra)
library(ggrepel)
library(GGally)
library(edgeR)
library(ineq)
library(dplyr)
library(stringr)
library(ComplexHeatmap)
library(MAGeCKFlute)



### Initialize ###

#directories
base_dir <- "./SSEA5_differentiation_screen/"
output_dir <- "./SSEA5_differentiation_screen/screen_analysis/"

#colormaps
ScatterPlotColorPanel = rev(c("#F8FA0D", "#F6DA23", "#F8BA43","#A5BE6A","#2DB7A3","#1389D2","#0262E0","#343DAE","#352A86"))

#functions
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}


### Pre-processing ###

#compile raw sgRNA counts for all samples 
sgrnas_lib <- fread("./epifactors_sgRNA_library/EpiFactors_CRISPRi_lib_v3.tsv", data.table = FALSE)
sgrna_counts <- sgrnas_lib[, c("sgRNA_ID", "sgRNA_Sequence", "Gene_Symbol")]
counts_folder <- paste0(base_dir, "sgRNA_counts")

for (val in list.files(path = counts_folder, pattern = ".rds")) {
  sample_name <- substr(val, 1, nchar(val)-11)
  sgrna_counts[,sample_name] <- readRDS(file = paste(counts_folder, val, sep = "/"))$sample
}

#CPM normalization of raw counts 
cpm.df <- cpm(sgrna_counts[,4:ncol(sgrna_counts)]+1) 
sgrna_cpms <- sgrna_counts
sgrna_cpms[,4:ncol(sgrna_cpms)] <- as.data.frame(cpm.df)

#save sgRNA counts and CPM matrices
saveRDS(sgrna_counts, file = paste0(base_dir, "sgRNA_counts_table.rds"))
saveRDS(sgrna_cpms, file =  paste0(base_dir, "sgRNA_cpm_table.rds"))

#output sgRNA counts and CPM matrices for supplemental tables 
write.table(sgrna_counts, file = paste0(base_dir, "sgRNA_counts_table.tsv"), 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

write.table(sgrna_cpms, file = paste0(base_dir, "sgRNA_cpm_table.tsv"), 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)



### Quality Control ###

#sgRNA abundance in plasmid versus D0 iPS gDNA libraries 
sgrna_cpms$density <- get_density(sgrna_cpms$plasmid_lib, sgrna_cpms$D0_iPS_1157.2_P38, n = 100)

pdf(paste0(output_dir, "qc_cpm_plasmid_v_D0_iPS_gDNA_libs.pdf"), width = 4, height = 4, useDingbats = FALSE)

gg <- ggplot(sgrna_cpms) + 
  geom_abline(slope = 1, intercept = 0, linetype = 2) + 
  geom_point(aes(log10(plasmid_lib), log10(D0_iPS_1157.2_P38), color = density)) + 
  theme_bw() + 
  scale_color_gradientn(colours = ScatterPlotColorPanel) +
  scale_x_continuous(limits = c(0, 3.5)) + 
  scale_y_continuous(limits = c(0, 3.5)) + 
  labs(x = "Plasmid Library log10(CPM)", y = "D0 iPS gDNA Library log10(CPM)") + 
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none")
ggMarginal(gg, type = "histogram", fill = "grey")

dev.off()

lib_cor <- cor(sgrna_cpms$plasmid_lib, sgrna_cpms$D0_iPS_1157.2_P38, method = "pearson")
print(lib_cor)
#r = 0.835346

#lorenz plot of the plasmid library 
neworder <- sgrna_cpms[order(sgrna_cpms$plasmid_lib),]

lcolc <- Lc(sgrna_cpms$plasmid_lib)
lcdf <- data.frame(L = rev(1-lcolc$L), p = lcolc$p, Uprob = c(1:length(lcolc$L)/length(lcolc$L)))

pdf(paste0(output_dir, "qc_lorenz_curve_plasmid_lib.pdf"), width = 4, height = 4, useDingbats = FALSE)

ggplot(lcdf, aes(x = Uprob, y = L)) + 
  geom_line(colour = hcl(h=15, l=65, c=100)) + 
  geom_line(aes(x = p, y = p)) +
  scale_x_continuous(breaks=c(0,0.25,0.5,0.75,1)) + scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1)) +
  theme_bw() + 
  labs(title="Lorenz Curve of sgRNA Representation\nPlasmid Library", x="Descending sgRNA Abundance Rank", y = "Cumulative Fraction of Reads") +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))

dev.off()

gini <- ineq(sgrna_cpms$plasmid_lib, type = "Gini")
auc <- (gini+1)*0.5
print(auc)
#auc = 0.6239564

qs <- quantile(sgrna_cpms$plasmid_lib, c(.1,.9))
skewness_ratio <- qs[2]/qs[1]
print(skewness_ratio)
#skew ratio = 3.481969

#lorenz plot of the D0 iPS gDNA library 
neworder <- sgrna_cpms[order(sgrna_cpms$D0_iPS_1157.2_P38),]

lcolc <- Lc(sgrna_cpms$D0_iPS_1157.2_P38)
lcdf <- data.frame(L = rev(1-lcolc$L), p = lcolc$p, Uprob = c(1:length(lcolc$L)/length(lcolc$L)))

pdf(paste0(output_dir, "qc_lorenz_curve_D0_iPS_gDNA_lib.pdf"), width = 4, height = 4, useDingbats = FALSE)

ggplot(lcdf, aes(x = Uprob, y = L)) + 
  geom_line(colour = hcl(h=15, l=65, c=100)) + 
  geom_line(aes(x = p, y = p)) +
  scale_x_continuous(breaks=c(0,0.25,0.5,0.75,1)) + scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1)) +
  theme_bw() + 
  labs(title="Lorenz Curve of sgRNA Representation\nD0 iPS Genomic DNA", x="Descending sgRNA Abundance Rank", y = "Cumulative Fraction of Reads") +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))

dev.off()

gini <- ineq(sgrna_cpms$D0_iPS_1157.2_P38, type = "Gini")
auc <- (gini+1)*0.5
print(auc)
#auc = 0.6060972

qs <- quantile(sgrna_cpms$D0_iPS_1157.2_P38, c(.1,.9))
skewness_ratio <- qs[2]/qs[1]
print(skewness_ratio)
#skew ratio = 2.691088

#plot pairwise sgRNA abundance across all samples
ggplot2::theme_set(ggplot2::theme_bw())

pdf(paste0(output_dir, "qc_pairwise_sgrna_abundance.pdf"), width = 10, height = 10, useDingbats = FALSE)
ggpairs(log10(sgrna_cpms[,4:(ncol(sgrna_cpms)-1)]))
dev.off()

#determine concordance of sgRNA abundances across screen replicates
res <- cor(sgrna_cpms[,c("D9_rep1_SSEA5_hi", "D9_rep1_SSEA5_low", "D9_rep2_SSEA5_hi", "D9_rep2_SSEA5_low", "D9_rep3_SSEA5_hi", "D9_rep3_SSEA5_low")])

column_ha = HeatmapAnnotation(Sample = factor( str_split_fixed(colnames(res), pattern = "_", n=3)[,3], levels = c("SSEA5_low","SSEA5_hi") ), 
                              border = TRUE)

hm <- Heatmap(res, #col=circlize::colorRamp2(c(0,1), c("white", "dark blue")),
              show_row_names = FALSE, 
              show_column_names = FALSE, 
              top_annotation = column_ha,
              show_row_dend = FALSE,
              border = TRUE, 
              column_title = "Concordance of Replicate Screens",
              rect_gp = gpar(col = "white", lwd = 1),
              name = "Pearson Correlation")

pdf(paste0(output_dir, "qc_D9_SSEA5_CRISPRi_screen_replicate_concordance.pdf"), width = 4.5, height = 3.5)  
draw(hm, merge_legend = TRUE)
dev.off()



### iPS SSEA-5 D9 CRISPRi Screen Analysis ###

#import MAGeCK results 
screen_name <- "iPS_SSEA5_D9_CRISPRi"
mageck_dir <- paste0(base_dir, "mageck_output")

file1 = file.path(paste0(mageck_dir, "/test/", screen_name, ".gene_summary.txt"))
file2 = file.path(paste0(mageck_dir, "/test/", screen_name, ".sgrna_summary.txt"))

grra <- ReadRRA(file1)
sgrra = ReadsgRRA(file2)

gene_summary <- fread(file1, data.table = FALSE)

#calculate sgRNA-level CRISPR scores
sgrna_cpms$CS_rep1 <- log2(sgrna_cpms$D9_rep1_SSEA5_low / sgrna_cpms$D9_rep1_SSEA5_hi)
sgrna_cpms$CS_rep2 <- log2(sgrna_cpms$D9_rep2_SSEA5_low / sgrna_cpms$D9_rep2_SSEA5_hi)
sgrna_cpms$CS_rep3 <- log2(sgrna_cpms$D9_rep3_SSEA5_low / sgrna_cpms$D9_rep3_SSEA5_hi)

#rename non-targeting control sgRNAs so that they aren't collapsed 
inds <- sgrna_cpms$Gene_Symbol == "NO-TARGET"
sgrna_cpms[inds,"Gene_Symbol"] <- paste0(sgrna_cpms[inds,"Gene_Symbol"], c(1:sum(inds)))

#calculate gene-level CRISPR scores 
df.cs <- data.frame(gene = sgrna_cpms$Gene_Symbol %>% unique())
df.cs$cs_rep1 <- aggregate(sgrna_cpms$CS_rep1, list(sgrna_cpms$Gene_Symbol), FUN=mean)$x
df.cs$cs_rep2 <- aggregate(sgrna_cpms$CS_rep2, list(sgrna_cpms$Gene_Symbol), FUN=mean)$x
df.cs$cs_rep3 <- aggregate(sgrna_cpms$CS_rep3, list(sgrna_cpms$Gene_Symbol), FUN=mean)$x
rownames(df.cs) <- df.cs$gene

#isolate hits from MAGeCK pipeline (log2 fold change > 1 and FDR < 0.05)
df.hits <- gene_summary %>% filter(`pos|fdr` < 0.05 & `pos|lfc` > 1)

#visualize MAGeCK-defined screen hits
df.cs$highlight <- "Not a hit"
df.cs[ df.cs$gene %in% df.hits$id, "highlight"] <- "Hit"
df.cs[grepl("NO-TARGET", df.cs$gene), "highlight"] <- "Non-targeting control"
df.cs$highlight <- factor(df.cs$highlight, levels = c("Not a hit", "Non-targeting control", "Hit"))
df.cs <- df.cs[is.finite(rowSums(df.cs[,c("cs_rep1", "cs_rep2", "cs_rep3")])),]

pdf(paste0(output_dir, "D9_SSEA5_CRISPRi_screen_log2FC_rep1_v_rep2.pdf"), width = 7, height = 6, useDingbats = FALSE)

ggplot(df.cs %>% arrange(highlight), aes(x = cs_rep1, y = cs_rep2, color = highlight)) + 
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  geom_vline(xintercept=0, linetype="dashed", color = "black") + 
  geom_point() + 
  scale_color_manual(values=c("grey", "red", "blue")) + 
  labs(x = "Replicate 1, log2(Fold Change)", y = "Replicate 2, log2(Fold Change)") + 
  #xlim(-4,4) + ylim(-3,3) + 
  ggtitle("D9 iPS SSEA-5 CRISPRi Screen") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title= element_blank()
  ) +
  geom_label_repel( 
    data=df.cs %>% filter(gene %in% df.hits$id), 
    aes(label=gene),
    #nudge_x = .1, nudge_y = .1,
    label.size = NA,
    fill = alpha(c("white"),0.1),
    box.padding   = 0.35, 
    point.padding = 0.5,
    max.overlaps = 100,
    segment.color = "black"
    #position = position_dodge(0.8)
  ) 

dev.off()

#STRING analysis of hits 
library(STRINGdb)

#full STRING network
string_db <- STRINGdb$new(version = "11.5", 
                          species = 9606, 
                          score_threshold = 400, 
                          network_type = "full", 
                          input_directory = "")
mapped <- string_db$map(df.hits, "id", removeUnmappedRows = TRUE)
string_hits <- mapped$STRING_id

pdf(paste0(output_dir, "D9_SSEA5_CRISPRi_screen_STRING_full_network_score400_all_hits.pdf"), width = 5, height = 5, useDingbats = FALSE)
string_db$plot_network(string_hits)
dev.off()

clustersList <- string_db$get_clusters(mapped$STRING_id)

pdf(paste0(output_dir, "D9_SSEA5_CRISPRi_screen_STRING_full_network_score400_clusters.pdf"), width = 5, height = 8, useDingbats = FALSE)
par(mfrow=c(3,2))
for(i in seq(1:5)){
  string_db$plot_network(clustersList[[i]])
}
dev.off()

enrichment <- string_db$get_enrichment( string_hits )
c(11,15,16,17,18,19,21,22,24,25,29)

#physical protein interaction network
string_db <- STRINGdb$new(version = "11.5", 
                          species = 9606, 
                          score_threshold = 400, 
                          network_type = "physical",
                          input_directory = "")
mapped <- string_db$map(df.hits, "id", removeUnmappedRows = TRUE)
string_hits <- mapped$STRING_id

pdf(paste0(output_dir, "D9_SSEA5_CRISPRi_screen_STRING_physical_network_score400_all_hits.pdf"), width = 5, height = 5, useDingbats = FALSE)
string_db$plot_network(string_hits)
dev.off()

clustersList <- string_db$get_clusters(mapped$STRING_id)

pdf(paste0(output_dir, "D9_SSEA5_CRISPRi_screen_STRING_physical_network_score400_clusters.pdf"), width = 5, height = 8, useDingbats = FALSE)
par(mfrow=c(3,2))
for(i in seq(1:6)){
  string_db$plot_network(clustersList[[i]])
}
dev.off()







#perform enrichment analysis on screen hits 
genelist = grra$Score
names(genelist) = grra$id
genelist = sort(genelist, decreasing = TRUE)
head(genelist)
hgtRes1 = EnrichAnalyzer(genelist, method = "HGT", 
                         type = "Pathway", organism = "hsa")
head(hgtRes1@result)

barplot(hgtRes1, showCategory = 5)

gseRes1 = EnrichAnalyzer(genelist, method = "GSEA", type = "Pathway", organism = "hsa")

idx = which(gseRes1$NES>0)[1]
gseaplot(gseRes1, geneSetID = idx, title = gseRes1$Description[idx])




#epifactors_main <- fread(file = "./epifactors_sgRNA_library/EpiGenes_main_v1.7.3.tsv", data.table = FALSE)
#rownames(epifactors_main) <- epifactors_main$HGNC_symbol

#df.epifunctions <- tidyr::separate_rows(epifactors_main[,c("Function", "HGNC_symbol")], Function, sep=',') %>% as.data.frame()
#df.epifunctions <- tidyr::separate_rows(epifactors_main[,c("Modification", "HGNC_symbol")], Modification, sep=',') %>% as.data.frame()
#df.epifunctions <- tidyr::separate_rows(epifactors_main[,c("Complex_name", "HGNC_symbol")], Complex_name, sep=',') %>% as.data.frame()
#df.epifunctions <- tidyr::separate_rows(epifactors_main[,c("Target", "HGNC_symbol")], Target, sep=',') %>% as.data.frame()
#df.epifunctions <- tidyr::separate_rows(epifactors_main[,c("Specific_target", "HGNC_symbol")], Specific_target, sep=',') %>% as.data.frame()
#colnames(df.epifunctions) <- c("term", "gene")

#ans.tf <- enricher(df.hits$id, TERM2GENE=df.epifunctions)






#plot sgRNA ranks for select genes 
ino80 <- c("ACTR5", "INO80C", "INO80E", "UCHL5", "NFRKB")
prc <- c("RNF2", "EZH2", "PCGF3", "ASXL2", "CBX2", "L3MBTL2")
heterochromatin <- c("EHMT2", "EHMT1", "SETDB1", "ATF7IP", "CBX3", "MBD6", "MPHOSPH8")
saga <- c("TAF6L", "TAF5L", "KAT2A", "TADA2B", "CCDC101", "SUPT7L", "TADA3", "TADA1")
#e2f6 <- c("EHMT2", "EZH2", "EHMT1", "CBX3", "RNF2", "L3MBTL2")

sgRankView(sgrra, 
           gene = c(ino80, prc, heterochromatin, saga), 
           top = 0, 
           bottom = 0, 
           neg_ctrl = "NO-TARGET", 
           bg.col = "white", 
           width = 2.75,
           height = 5, 
           filename = paste0(output_dir, "sgRNA_rank_view_hits.pdf")
)




#enrichment analysis 
geneList <- grra$Score
names(geneList) <- grra$id
enrich = EnrichAnalyzer(geneList = geneList[geneList>0.5], 
                        method = "HGT", 
                        type = "KEGG+REACTOME+GOBP+GOMF+GOCC", 
                        filter = TRUE)

EnrichedView(enrich, mode = 1, top = 20)
EnrichedView(enrich, mode = 2, top = 20)

enrich = EnrichAnalyzer(geneList = geneList[geneList>0.5], 
                        method = "HGT", 
                        type = "Pathway", 
                        organism = "hsa")





enrich = EnrichAnalyzer(geneList = geneList[geneList>0.5], 
                        method = "ORT", 
                        type = "KEGG+REACTOME+GOBP+GOMF+GOCC",
                        filter = TRUE)






library(clusterProfiler)
tmp <- gene_summary %>% filter(`pos|fdr` < 0.05 & `pos|lfc` > 0.5)


ego <- enrichGO(gene          = tmp$id,
                #universe      = df.mageck$id,
                OrgDb         = org.Hs.eg.db,
                keyType       = "SYMBOL",
                ont           = "ALL", #BP, MF, CC
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
)
head(ego)

cnetplot(ego,4)
barplot(ego)


geneList = gene_summary$`pos|lfc`
names(geneList) = as.character(gene_summary$id)
geneList = sort(geneList, decreasing = TRUE)

ego3 <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              ont          = "BP",
              keyType      = "SYMBOL",
              #minGSSize    = 100,
              #maxGSSize    = 500,
              #pvalueCutoff = 0.05,
              verbose      = FALSE)


gene <- names(geneList)[geneList > 1]
df.bitr <- bitr(gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

kk <- enrichKEGG(gene         = df.bitr$ENTREZID,
                 organism     = 'hsa',
                 keyType      = "kegg",
                 pvalueCutoff = 0.05)
head(kk)

browseKEGG(kk, 'hsa04110')

library(ReactomePA)
x <- enrichPathway(gene         = df.bitr$ENTREZID, 
                   pvalueCutoff = 0.05, 
                   readable     = TRUE)

head(x)













