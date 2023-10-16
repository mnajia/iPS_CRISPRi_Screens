#Project: iPS CRISPRi
#plot_epifactors_library_composition.R
#Objective: Visualize the composition of genes in the EpiFactors CRISPRi library
#Mohamad Najia

library(ggpubr)
library(ggplot2)
library(stringr)
library(dplyr)
library(tidyverse)
library(data.table)


### Initialization ###

#Import the library annotation file 
epifactors_main <- fread(file = "./epifactors_sgRNA_library/EpiGenes_main_v1.7.3.tsv", data.table = FALSE)


### Visualize a summary of gene functions within the EpiFactors library ###

#quantify epigenetic functions of genes within the library 
df.functions <- data.frame(Function = unlist(str_split(epifactors_main$Function, ", ")))
df.function.counts <- aggregate(df.functions, list(Term = df.functions$Function), length)
df.function.counts <- df.function.counts[-c(14),]

#plot a summary of gene functions in the Epifactors Library
pdf("./epifactors_sgRNA_library/epifactors_crispri_library_epigenetic_functions_distribution.pdf", width = 6, height = 5, useDingbats = FALSE)
par(cex.main=0.8,mar=c(2,2,2,2))

ggdotchart(df.function.counts, x = "Term", y = "Function",
           color = "#FC4E07", 
           sorting = "descending",
           add = "segments",
           ylab = "Number of Genes in Library",
           xlab = "Function",
           rotate = TRUE,
           dot.size = 7, 
           label = round(df.function.counts$Function),
           font.label = list(color = "white", size = 9, 
                             vjust = 0.5), 
           ggtheme = theme_pubr() 
)

dev.off()


### Visualize chromatin complexes represented within the library ###

#import chromatin complex annotation file and create dataset 
df.complexes <- fread(file = "./epifactors_sgRNA_library/EpiGenes_complexes_v1.7.3.tsv", data.table = FALSE)
complex.members <- str_split(df.complexes$UniProt_ID, ",")
df.complexes$counts <- as.numeric(lapply(complex.members, length))
data <- data.frame(individual = df.complexes$Complex_name, 
                   group = as.factor(df.complexes$Group_name), 
                   value = df.complexes$counts)
data <- data %>% arrange(group, value)

empty_bar <- 3
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$group), ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$group <- rep(levels(data$group), each=empty_bar)
data <- rbind(data, to_add)
data <- data %>% arrange(group)
data$id <- seq(1, nrow(data))

label_data <- data
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

base_data <- data %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]


#plot a summary of epigenetic complexes represented within the library
pdf("./epifactors_sgRNA_library/epifactors_crispri_library_epigenetic_complexes_distribution.pdf", width = 8, height = 8, useDingbats = FALSE)

p <- ggplot(data, aes(x=as.factor(id), y=value, fill=group)) + 
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
  
  geom_segment(data=grid_data, aes(x=end, y=40, xend=start, yend=40), colour="grey", alpha=1, size=0.3, inherit.aes=FALSE) +
  geom_segment(data=grid_data, aes(x=end, y=30, xend=start, yend=30), colour="grey", alpha=1, size=0.3, inherit.aes=FALSE) +
  geom_segment(data=grid_data, aes(x=end, y=20, xend=start, yend=20), colour="grey", alpha=1, size=0.3, inherit.aes=FALSE) +
  geom_segment(data=grid_data, aes(x=end, y=10, xend=start, yend=10), colour="grey", alpha=1, size=0.3, inherit.aes=FALSE) +
  
  annotate("text", x=rep(max(data$id),4), y=c(10, 20, 30, 40), label=c("10","20","30","40") , color="grey", size=2, angle=0, fontface="bold", hjust=1) +
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
  ylim(-100,120) +
  theme_minimal() +
  theme(
    legend.position="none",
    axis.text=element_blank(),
    axis.title=element_blank(),
    panel.grid=element_blank(),
    plot.margin=unit(rep(-1,4), "cm") 
  ) +
  coord_polar() + 
  geom_text(data=label_data, aes(x=id, y=value+10, label=individual, hjust=hjust), color="black", fontface="bold", alpha=0.6, size=2.5, angle=label_data$angle, inherit.aes=FALSE) +
  geom_segment(data=base_data, aes(x=start, y=-5, xend=end, yend=-5), colour="black", alpha=0.8, size=0.6, inherit.aes=FALSE) +
  geom_text(data=base_data, aes(x=title, y=-18, label=group), colour="black", alpha=0.8, size=2, fontface="bold", inherit.aes=FALSE)

p

dev.off()




