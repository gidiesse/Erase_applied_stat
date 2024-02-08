# genome_eda: in this script we began to explore the PCA scores. In particular 
# we explored the relations between the first three scores, also differentiating 
# between modern and ancient samples and by main ancestry. 

# Libraries and wd 
library(bigsnpr) # for loading bed bim fam files
library(ggplot2) 

# Set directory with the path to the folder of the github repo
setwd("~/Desktop/ERASE_project2023_desgrp")

# If it's the first time you read the .bed file, run these two lines
bed_file <- "/data/snp/adna_pca_newiid_2convertV3.bed"
rds <- snp_readBed(bedfile)

# If you already read the .bed file once, run these line instead
rds <- "data/snp/adna_pca_newiid_2convertV3.rds"

# Loading the data from backing files (snp object loaded in the environment)
genetic <- snp_attach(rds)
big_counts(genetic$genotypes, ind.col = 1:10)

### PCA qualitative analysis ###

# Build a header vector with names for the columns of PC scores file
col_names = c("sample")
for (n in 1:100){
  col_names = c(col_names, sprintf("score_%s", n))
}
col_names = c(col_names, "pop")

# Load the PC scores
pc_file <- "data/PCA/adna_pca_newiid_2convertV3_100.evec"
evecDat <- read.table(pc_file, col.names=col_names, fill = TRUE)

# Load the variability correspondent to each PC
var_file <- "data/PCA/adna_pca_newiid_2convertV3_100.eval"
var.values <- read.table(var_file)$V1

# Scree plot

var.explained <- var.values / sum(var.values)
pc.num <- seq(1, length(var.explained))
scree <- data.frame(pc.num, var.explained)
quartz()
ggplot(data=scree, aes(pc.num,var.explained)) + geom_line() + geom_point() +
  theme(panel.background = element_rect(fill = "#cae8ff",
                                        colour = "#cae8ff",
                                        size = 0.5, linetype = "solid"),
        plot.title = element_text(hjust = 0.5)) +
  labs(y="Variance explained", x="Number of PCs") +
  geom_vline(xintercept = 100, linetype='longdash', color='#ff0000') +
  ggtitle("Screeplot")

cumulative.var.explained <- cumsum(var.values) / sum(var.values)
elbow.plot <- data.frame(pc.num, cumulative.var.explained)
quartz()
ggplot(data=elbow.plot, aes(pc.num,cumulative.var.explained)) + geom_line() + geom_point()

# Load the metadata
metadata <- read.csv("data/erase_metadata.csv", header=TRUE)
load("~/Downloads/MA_coord_df.Rda")
metadata$Macroarea = points_three$Macroarea
metadata$lat_num = points_three$lat_num
metadata$long_num = points_three$long_num
sum(is.na(metadata$date_fix))
# Join the two dataframes to get an aggregate: we can now use the metadata 
# information to help interpretation of PC scores 
pca.metadata <- merge(evecDat, metadata, by="sample") # unique identifier is sample

# plots

# (1): score2 vs score1, colors on ancient/modern individual. 
quartz()
ggplot(data=pca.metadata, aes(score_1, score_2, colour=project)) + 
  geom_point(size=0.5)

# (2): score3 vs score1, colors on ancient/modern individual
quartz()
ggplot(data=pca.metadata, aes(score_1, score_3, colour=project)) + 
  geom_point(size=0.5)

# (3): score3 vs score1, transformed in polar coordinates, colors on ancient/modern individual
quartz()
centered_score_1 = pca.metadata$score_1 - mean(pca.metadata$score_1)
centered_score_2 = pca.metadata$score_2 - mean(pca.metadata$score_2)
pca.metadata.plot = data.frame(centered_score_1, centered_score_2, pca.metadata$project)

ggplot(data=pca.metadata.plot, aes(sqrt(centered_score_1^2+centered_score_2^2), atan2(centered_score_2, centered_score_1), colour=project)) + 
  geom_point(size=0.5)

ggplot(data=pca.metadata, aes(sqrt(score_1^2+score_3^2), atan2(score_3, score_1), colour=project)) + 
  geom_point(size=0.5)

# (4): score2 vs score1, colors on Macroarea individual, NAs in grey
quartz()
ggplot(data=pca.metadata, aes(score_1, score_2, colour=Macroarea)) + 
  geom_point(size=0.5)

# (5): score2 vs score 1, colors on Macroarea, NAs not considered for clarity
quartz()
ggplot(data=subset(pca.metadata, !is.na(Macroarea)), aes(score_2, score_1)) + 
  geom_point(size=0.5, aes(color=Macroarea)) 

# (6):score3 vs score 1, colors on Macroarea, NAs in grey
quartz()
ggplot(data=pca.metadata, aes(score_1, score_3, colour=Macroarea)) + 
  geom_point(size=0.5)

# (7): score3 vs score 1, colors on Macroarea, NAs not considered for clarity
quartz()
ggplot(data=subset(pca.metadata, !is.na(Macroarea)), aes(score_1, score_3, colour=Macroarea)) + 
  geom_point(size=0.5)

# We now divide our data in two samples: modern sample and ancient sample
pca.metadata.ancient <- subset(pca.metadata, project=="ancient")
not.anc <- c("Iran_Neolithic", "Levant_N", "out", "Steppe")
idx.in <- (pca.metadata.ancient$Main_Ancestry!=not.anc[1] &
              pca.metadata.ancient$Main_Ancestry!=not.anc[2] &
              pca.metadata.ancient$Main_Ancestry!=not.anc[3] &
              pca.metadata.ancient$Main_Ancestry!=not.anc[4]) | is.na(pca.metadata.ancient$Main_Ancestry)
pca.metadata.ancient <- pca.metadata.ancient[idx.in,]
pca.metadata.ancient$Main_Ancestry <- factor(pca.metadata.ancient$Main_Ancestry,
                                             levels=c("WHG", "EHG", "Anatolia_N", 
                                                      "CHG", "North_Africa"))
# (8): score2 vs score1, colors on Main ancestry for 
quartz()
ggplot(data=subset(pca.metadata.ancient, !is.na(Main_Ancestry)), aes(score_1, score_2, colour=Main_Ancestry)) + geom_point(size=0.5)
quartz()
ggplot(data=pca.metadata.ancient[is.na(pca.metadata.ancient$Main_Ancestry),], aes(score_1, score_2, colour=Main_Ancestry)) + 
  geom_point(size=0.5) +
  geom_point(data=pca.metadata.ancient[!is.na(pca.metadata.ancient$Main_Ancestry),],  aes(score_1, score_2, colour=Main_Ancestry)) +
  scale_color_manual(values = c('red', 'darkorange', 'yellow', 'green', 'mediumblue')) +
  labs(y="Second Score", x="First Score") + 
  theme( panel.background = element_rect(fill = "#cae8ff",
                                       colour = "#cae8ff",
                                       size = 0.5, linetype = "solid"))
  
  
dim(pca.metadata.ancient[!is.na(pca.metadata.ancient$Main_Ancestry),])


# 3rd plot
check <- subset(pca.metadata, project=="modern")
# not feasible to plot wrt pop, too many categories
quartz()
ggplot(data=subset(check, !is.na(pop.x)), aes(score_1, score_2, colour=pop.x)) + geom_point(size=0.5)



