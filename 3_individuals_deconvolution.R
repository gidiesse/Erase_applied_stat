# In this file we explored different ways of assigning ancestries to individuals who didn't have one 
# (6513 of the 6580 individuals in our dataset do not have an ancestry assigned)
# We tried three different approaches. The first 2 were Gaussian Mixture Models and a C-means clustering
# (having such few samples with labeled ancestries - only 67 out of 6580 - we needed to use 
# unsupervised methods). Both of these had disastrous results and so were quickly abandoned. 
# Our last idea was to do a sort of deconvolution. Essentially we divided the 67 individuals into their 
# own respective ancestries. We then computed the mean of the pca scores for each ancestry resulting 
# into 5 individuals that were the "representatives" of their ancestry. These would be our centroids. 
# After that, we took the remaining 6513 individuals without an ancestry and computed the cosine similarity
# between their PCA scores and the PCA scores of each of the 5 "centroid" individuals. We then used a 
# softmax to normalize the similarity scores so that now each individual had 5 scores all between 0 and 1
# which summed to 1 - essentially each score represented how much the individual belonged to each ancestry. 
# We then assigned each individual to the ancestry that had the biggest score. This approach worked rather 
# well. In fact of the 67 individuals that already had ancestry labels, using this method 66 were assigned 
# corrrectly. 


### COMPOSITION OF INDIVIDUALS - prototype ###
# For each individual who belong to the Iron Age and Bronze Age period, 
# we want to obtain 5 numbers that add up to 1 (probabilties in a loose sense),
# which represent his belonging to each of the 5 main ancestries we considered
# ANCESTRIES CONSIDERED -> (CHG, EHG, WHG, Anatolia_n, North_Africa)

# Input:  representation of each individual in the coordinates of the first 100 PCs
# Output: composition of the individuals, based on the Macroancestry. A higher
#         score means a higher similarity to the centroid of that particular Macroancestry. 

# for cosine similarity
library(lsa)

# for c-means
library(e1071)

# for Gaussian mixture models
library(ClusterR)

# for data wrangling
library(dplyr)

# for normality checks
library(mvtnorm)
library(MVN)

# Set directory with the path to the folder of the github repo
setwd("~/Desktop/ERASE_project2023_desgrp")

# Build a header vector with names for the columns of PC scores file
col_names = c("sample")
for (n in 1:100){
  col_names = c(col_names, sprintf("score_%s", n))
}
col_names = c(col_names, "pop")

# Load the PC scores
pc_file <- "data/PCA/adna_pca_newiid_2convertV3_100.evec"
evecDat <- read.table(pc_file, col.names=col_names, fill=TRUE)

# Load the metadata
metadata <- read.csv("data/erase_metadata_interpolated.csv", header=TRUE)

# Join the two dataframes to get an aggregate: we can now use the metadata 
# information to help interpretation of PC scores 
pca.metadata <- merge(evecDat, metadata, by="sample") # unique identifier is sample

# remove the modern individuals (we agreed on focusing only on ancient ones)
pca.metadata.old <- pca.metadata[c(which(pca.metadata$project=='ancient')),]

# we check how many individuals we have for each ancestry
ancestries <- c() 
for(anc in unique(pca.metadata.old$Main_Ancestry)){
  ancestries <- c(ancestries, nrow(pca.metadata.old[which(pca.metadata.old$Main_Ancestry==anc),]))
}
ancestries 

# Since we have more than one individual for each main Ancestry that we want to consider,
# we decided to compute the centroids as the average of the individuals belonging
# to that ancestry
our_ancestries <- c('WHG', 'EHG', 'Anatolia_N', 'CHG', 'North_Africa')
centroids <- c()
center.individuals <- c()
#ancestries.center <- c()
for(anc in our_ancestries){
  individuals <- pca.metadata.old[which(pca.metadata.old$Main_Ancestry==anc),]
  centroids <- cbind(centroids, sapply(individuals[,2:101], mean))
  #center.individuals <- rbind(center.individuals, individuals)
  center.individuals <- c(center.individuals, which(pca.metadata.old$Main_Ancestry==anc))
  #ancestries.center <- rbind(ancestries.center, rep(anc, nrow(individuals)))
}
colnames(centroids) <- our_ancestries
centroids <- data.frame(centroids)
# In centroids we have 5 columns, that correspond to the centroid of each ancestry
# (So the columns are 100 entries vector)

# Keeping observations in Iron Age and Bronze Age and filter out the individuals
# we used to compute the centroids
#time_filter <- c('BronzeAge', 'IronAge')
#pca.metadata.final <- pca.metadata.old[which(pca.metadata.old$Macro_Period==time_filter),]
# we decide to keep all ancient individuals
#pca.metadata.final <- pca.metadata.old[which(is.na(pca.metadata.old$Main_Ancestry)),]
#pca.metadata.final <- rbind(pca.metadata.final, pca.metadata.old[center.individuals,])
pca.metadata.final <- pca.metadata.old
# now compute the cosine similarity between each observation and the centroids
sim.matrix <- c()
for(anc in our_ancestries){
  similarity <- c()
  for(idx in 1:nrow(pca.metadata.final)){
      similarity <- c(similarity, cosine(as.numeric(pca.metadata.final[idx, 2:101]), centroids[,anc]))
  }
  sim.matrix <- cbind(sim.matrix, similarity)
}
colnames(sim.matrix) <- our_ancestries
# in sim.matrix we have for each observation 5 numbers included in [-1,1] that 
# correspond to the similarity of that individual to each of the 5 centroids
# The closer the cosine similarity is to 1, the closer the individual is to the centroid
# If the cosine similarity is equal to zero, the vectors are orthogonal
# If the cosine similarity is -1, the vectors are one the opposite of the other.

# We use softmax function to normalize the similarity scores and obtaining 
# the "probabilities" of the composition of each individual.
softmax <- function(x) {
  norm_constant <- sum(exp(x))
  exp(x)/norm_constant
}

# Composition of the individuals: each row represents an individual and each column is an ancestry, each 
# value(i,j) of the matrix shows how much individual i belongs to ancestry j. Obviously the sum of each 
# row is 1.
comp.matrix <- c()
for(idx in 1:nrow(sim.matrix)){
  comp.matrix <- rbind(comp.matrix, softmax(sim.matrix[idx,]))
}
comp.matrix
comp.matrix <- as.data.frame(comp.matrix)

# Diagnostic of the composition matrix obtained: 
# We're happy if the Ancestry of the individuals is clearly distinguishable in terms of
# the composition we have obtained. 

# From the boxplot we see that there is certainly variability between the 
# composition of the ancestries, good start!
boxplot(comp.matrix, col=rainbow(5))

# We are in a semi-supervised setting, so let's check if the individuals for 
# which we already know the main ancestry (through biological tests) have been
# sorted correctly 

# We extract the individuals belonging to the known ancestries and check if they 
# have been assigned correctly
WHG.idx <- which(pca.metadata.final$Main_Ancestry=='WHG')
EHG.idx <- which(pca.metadata.final$Main_Ancestry=='EHG')
AnatoliaN.idx <- which(pca.metadata.final$Main_Ancestry=='Anatolia_N')
CHG.idx <- which(pca.metadata.final$Main_Ancestry=='CHG')
North_Africa.idx <- which(pca.metadata.final$Main_Ancestry=='North_Africa')
idxs <- list(WHG.idx, EHG.idx, AnatoliaN.idx, CHG.idx, North_Africa.idx)
assigned.ancestries <- matrix(0, 2, length(our_ancestries))
for (idx in seq(length(idxs))){
  assigned.anc <- apply(comp.matrix[idxs[[idx]],], 1, which.max)
  assigned.ancestries[1, idx] <- sum(assigned.anc==idx)
  assigned.ancestries[2, idx] <-  length(assigned.anc)
}
colnames(assigned.ancestries) <- our_ancestries
rownames(assigned.ancestries) <- c("Correctly assigned", "Total individuals")
assigned.ancestries
# NB the correspondence between 1-5 and the ancestries is as follows: 
# 1 - WHG, 2-EHG, 3-Anatolia_N, 4-CHG, 5-North_Africa

# We see that out of 67 individuals for which we had labeled ancestries, 66 have been assigned correctly

# to check the composition of the ancestry individuals "by hand"
comp.matrix[WHG.idx,]

# check composition of Steppe individuals, we expect high percentage of
# EHG and WHG and very little of North_Africa from our biological knowledge (Dr.Raveane)
unique(pca.metadata.final$Main_Ancestry)
idx.steppe <- which(pca.metadata.final$Main_Ancestry=='Steppe')
comp.matrix[idx.steppe,]
idx.iran <- which(pca.metadata.final$Main_Ancestry=="Iran_Neolithic")
comp.matrix[idx.iran,]

### We are very happy with these results, as it works almost perfectly. 
### Out of 67 ancestry individuals, the composition we obtained allowed us to 
### correctly classify 66 of them. This is an index of having obtained a good 
### deconvolution for our individuals.

### Let's look at how much dispersion we have for each ancestry composition.
### We would like to have a high dispersion, so that our individuals are 
### distinguishable based on the ancestry deconvolution.

# mean for the composition of each ancestry
m<-colMeans(comp.matrix)
m
# We see that WHG and EHG are the most prominent groups and CHG and North_africa
# are the least prominent. 

# standard deviation for the composition of each ancestry
s<- sqrt(sapply(comp.matrix, var))
s
# Let's compute confidence intervals, for composition of each ancient individuals 
#Since we have a large sample size we can use property of asymptotic Normality
#quantile di 0.975 = 1.96 
IC<- rbind(ifelse(m-1.96*s>0, m-1.96*s, 0), m,m+1.96*s)
rownames(IC) <- c('inf','mean','sup')
IC

plot(seq(1:5),m,ylim=c(0,0.5))
points(1:5, IC[2,], pch=16, col=1:5)
for(i in 1:5) segments(i, IC[1,i], i, IC[3,i], lwd=2, col=i)
points(1:5, IC[1,], pch='-', col=1:5)
points(1:5, IC[3,], pch='-', col=1:5)
title("Confidence Intervals for each Ancestry")

# To formally check if the means of the different ancestries are different, we 
# will run a MANOVA test where the response variable is the vector of 
# the compositions, while we obtain the groups (one-way Manova) by assigning 
# each individual to its most probable ancestry (by choosing the highest 
# composition score)

# STEP 1 --> assign each individual to its ancestry 
assigned_ancestry <- c()
for (i in 1:dim(comp.matrix)[1]) {
  assigned_ancestry <- c(assigned_ancestry, which.max(comp.matrix[i,]))
}
assigned_ancestry <- as.factor(assigned_ancestry)
comp.matrix.manova <- cbind(comp.matrix,assigned_ancestry)
count(comp.matrix.manova, assigned_ancestry)

# STEP 2 --> we check assumption of gaussianity
mvn(comp.matrix[assigned_ancestry == 1,])
mvn(comp.matrix[assigned_ancestry == 2,])
mvn(comp.matrix[assigned_ancestry == 3,])
mvn(comp.matrix[assigned_ancestry == 4,])
mvn(comp.matrix[assigned_ancestry == 5,])

# The assumption of gaussianity is not respected for any of the groups. We doubt 
# that the covariance matrices are the same. We decide nonetheless to carry on
# with the test. 

# STEP 3 --> MANOVA TIME <3
fit <- manova(as.matrix(comp.matrix) ~ assigned_ancestry)
summary.manova(fit, tol=0)
# we have had to set tol to 0 because the residuals have df = 4 < 5
# this implies collinearity, so we need to look into this 

cor(comp.matrix)
image(cor(comp.matrix))

pca.metada.composition <- cbind(pca.metadata.final[,102:116], comp.matrix.manova)
write.csv(pca.metada.composition, "./corr_convolution.csv")


#### C-means ####
#center.ancestries <- center.individuals$Main_Ancestry
#center.individuals <- as.matrix(center.individuals[,2:101])
center.individuals <- which(!is.na(pca.metadata.final$Main_Ancestry))
center.ancestries <- pca.metadata.final$Main_Ancestry[center.individuals]

individuals <- as.matrix(pca.metadata.final[,2:101])
#all.individuals <- rbind(center.individuals, individuals)


#centers <- t(as.matrix(centroids))[,1:100]
centers <- center.individuals[c(1,37, 40, 61, 63),]
set.seed(4)
cl <- cmeans(x=individuals, centers=5, iter.max=150,
             dist="euclidean", m=1.5)
# m=20, centers=centers, WORKS but just 1 iteration
composition_cmeans <- cl$membership
colnames(composition_cmeans) <- c("Cluster1", "Cluster2", "Cluster3", "Cluster4", "Cluster5")
labels_cmeans <- as.matrix(cl$cluster)
colnames(labels_cmeans) <- c("ClusterAssigned")
c_means_individuals <- cbind(pca.metadata.final, composition_cmeans, labels_cmeans)
c_means_individuals <- c_means_individuals[, 102:122]
write.csv(c_means_individuals, "./c_means_composition.csv")

mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

mode(cl$cluster[which(pca.metadata.final$Main_Ancestry=='EHG')])
mode(cl$cluster[which(pca.metadata.final$Main_Ancestry=='WHG')])
mode(cl$cluster[which(pca.metadata.final$Main_Ancestry=='Anatolia_N')])
mode(cl$cluster[which(pca.metadata.final$Main_Ancestry=='CHG')])
mode(cl$cluster[which(pca.metadata.final$Main_Ancestry=='North_Africa')])

# Gaussian mixture models
gmm <- GMM(individuals, 4, "eucl_dist", seed=17)
pr  <- predict(gmm, newdata=individuals)
composition_gmm <- gmm$Log_likelihood
colnames(composition_gmm) <- c("Cluster1", "Cluster2", "Cluster3", "Cluster4", "Cluster5")
labels_gmm <- c()
for(i in 1:1414){
  labels_gmm <- rbind(labels_gmm, which(gmm$Log_likelihood[i,]==max(gmm$Log_likelihood[i,])))
}
colnames(labels_gmm) = c("ClusterAssigned")
gmm_individuals <- cbind(pca.metadata.final, composition_gmm, labels_gmm)
gmm_individuals <- gmm_individuals[, 102:122]
write.csv(gmm_individuals, "./gmm_composition.csv")

mode(labels_gmm[which(pca.metadata.final$Main_Ancestry=='WHG')])
mode(labels_gmm[which(pca.metadata.final$Main_Ancestry=='EHG')])
mode(labels_gmm[which(pca.metadata.final$Main_Ancestry=='Anatolia_N')])
mode(labels_gmm[which(pca.metadata.final$Main_Ancestry=='CHG')])
mode(labels_gmm[which(pca.metadata.final$Main_Ancestry=='North_Africa')])

# standardized variables
compute_posterior <- function(log_likelihoods, weights) {
  (exp(log_likelihoods)*weights) / as.numeric((exp(log_likelihoods)%*%weights))
}
std.individuals <- center_scale(individuals, mean_center = T, sd_scale = T)
gmm <- GMM(std.individuals, 5, dist_mode = "maha_dist", seed_mode = "random_subset", 
           km_iter = 10, em_iter = 10)
pr  <- predict(gmm, newdata=std.individuals)
composition_gmm <- gmm$Log_likelihood
colnames(composition_gmm) <- c("Cluster1", "Cluster2", "Cluster3", "Cluster4", "Cluster5")
labels_gmm <- c()
for(i in 1:1414){
  composition_gmm[i,] <- compute_posterior(composition_gmm[i,], gmm$weights)
  labels_gmm <- rbind(labels_gmm, which.max(composition_gmm[i,]))
}
colnames(labels_gmm) = c("ClusterAssigned")
gmm_individuals <- cbind(pca.metadata.final, composition_gmm, labels_gmm)
gmm_individuals <- gmm_individuals[, 102:122]
#write.csv(gmm_individuals, "./gmm_composition.csv")



mode(labels_gmm[which(pca.metadata.final$Main_Ancestry=='WHG')])
mode(labels_gmm[which(pca.metadata.final$Main_Ancestry=='EHG')])
mode(labels_gmm[which(pca.metadata.final$Main_Ancestry=='Anatolia_N')])
mode(labels_gmm[which(pca.metadata.final$Main_Ancestry=='CHG')])
mode(labels_gmm[which(pca.metadata.final$Main_Ancestry=='North_Africa')])







