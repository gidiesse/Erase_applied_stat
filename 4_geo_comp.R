# This is the file where we use spatial statistics.  The crux of our project is understanding 
# how main ancestry evolved in the different macroareas over time. 
# The challenge was that during the deconvolution (previous R script) we generated compositional data,
# by adding the constraint that the values had to sum to 1 we lost a degree of freedom and as such are
# no longer in a Euclidean space and couldn't use any of the statistical methods we saw in class as
# in class we always oeprated under the assumption of being in a Euclidean space. 
# As a result, we performed an ILR transformation of our compositional data and that brought us back 
# into a Euclidean space. We then performed another PCA to find a 4-dimensional orthonormal basis. 
# We then performed our geospatial analysis and finally used the inverse ILR transformation to go 
# back to our starting space where we could then interpret the results. 

## Load spatial packages
library(sp)           ## Data management
library(lattice)      ## Data management
library(gstat)        ## Geostatistics
library(compositions) ## Compositional data management
library(cowplot)      ## To work with ggpplot
library(ggplot2)      ## To make nice graphs

# Set directory with the path to the folder of the github repo
setwd("~/Desktop/ERASE_project2023_desgrp")

# import data
individuals.deconvolution <- read.csv("corr_convolution.csv")

plot(individuals.deconvolution$PC2, individuals.deconvolution$PC1, col=as.factor(individuals.deconvolution$assigned_ancestry))
legend('topright', levels(as.factor(individuals.deconvolution$assigned_ancestry)),  fill=levels(as.factor(individuals.deconvolution$assigned_ancestry)),bty='n')
title ("Score1 vs Score2 by ancestry")

plot(individuals.deconvolution$PC2, individuals.deconvolution$PC1, col=as.factor(individuals.deconvolution$Macroarea))
#legend('topright', levels(as.factor(individuals.deconvolution$Macroarea)),  fill=levels(as.factor(individuals.deconvolution$Macroarea)),bty='n')

legend('topright', levels(as.factor(individuals.deconvolution$Macroarea)),  fill=levels(as.factor(individuals.deconvolution$Macroarea)),bty='n')


# drop NAs
na.lat_lon <- is.na(individuals.deconvolution$lat_num)
individuals.deconvolution <- individuals.deconvolution[!na.lat_lon,]
na.macroarea <- is.na(individuals.deconvolution$Macroarea)
individuals.deconvolution <- individuals.deconvolution[!na.macroarea,]
na.macroperiod <- is.na(individuals.deconvolution$Macro_Period)
individuals.deconvolution <- individuals.deconvolution[!na.macroperiod,]

set.seed(42)
n <- dim(individuals.deconvolution)[1]
individuals.deconvolution$lat_num <- individuals.deconvolution$lat_num + rnorm(n, 0, 0.025)
individuals.deconvolution$long_num <- individuals.deconvolution$long_num + rnorm(n, 0, 0.025)

# here we transform into a compositional data data type
x <- acomp(individuals.deconvolution[,17:21])

# here we calculate the pcs --> we get 4 and this makes sense because we go from Rd -> 
# Rd-1
pcx = princomp(x)
summary(pcx)

pcx$loadings
par(mfrow=c(2,2))
barplot(pcx$loadings[,1])
barplot(pcx$loadings[,2])
barplot(pcx$loadings[,3])
barplot(pcx$loadings[,4])


# Now we create a dataset with all the pcs + long, lat, macroarea and macroperiod
dat <- data.frame(pcx$scores[,1:4], Macro_Period=individuals.deconvolution$Macro_Period, 
             Macroarea = individuals.deconvolution$Macroarea)

# create the spatial data frame
proj.info <- "+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
coords <- cbind(individuals.deconvolution$long_num, individuals.deconvolution$lat_num)
individuals.spatial <- SpatialPointsDataFrame(
  coords=as.matrix(coords),
  data=dat, 
  proj4string = CRS(as.character(proj.info)))


# For instance: multivariate variogram
g = gstat(NULL, "Comp1", Comp.1 ~ Macroarea + Macro_Period, individuals.spatial)
g = gstat(g, "Comp2", Comp.2 ~ Macroarea + Macro_Period, individuals.spatial)
g = gstat(g, "Comp3", Comp.3 ~ Macroarea + Macro_Period, individuals.spatial)
g = gstat(g, "Comp4", Comp.4 ~ Macroarea + Macro_Period, individuals.spatial)
g 

# Variogram modeling
# Empirical estimate
vm = variogram(g, cutoff=1500)
plot(vm)

# Fit valid model
vm.fit = gstat::fit.lmc(vm, g, vgm(1, "Exp", 800, 1))
plot(vm, vm.fit)

# plot with ggplot:
df.variogram <- data.frame(dist=vm$dist[1:15], gamma=vm$gamma[1:15])

ggplot(df.variogram, aes(x=dist, y=gamma)) +
  geom_point(size=3.5, colour='#3a9bdb') +
  ylim(0, 0.0045) +
  theme(plot.title = element_text(hjust = 0.5), 
        panel.background = element_rect(fill = "#cae8ff",
                                        colour = "#cae8ff",
                                        size = 0.5, linetype = "solid")) +
  labs(y="Gamma", x="Geodetic distance") +
  ggtitle("Cross-Variogram of PC1 and PC4 in ILR space")

# Cross-variograms are not pure nuggets, so technically the four pcs are not 
# uncorrelated, but we will assume that they are anyways and fit a specific model
# for each pc

# Observations for checking the composition coefficients
mp <- unique(individuals.spatial$Macro_Period)
ma <- unique(individuals.spatial$Macroarea)
all.comb <- as.data.frame(expand.grid(Macroarea=ma, Macro_Period=mp))
coords.ac <- cbind(individuals.deconvolution$long_num[1:84], 
                   individuals.deconvolution$lat_num[1:84])
all.comb <- SpatialPointsDataFrame(
  coords=as.matrix(coords.ac),
  data=all.comb, 
  proj4string = CRS(as.character(proj.info)))

##### Here we focus on the model for PC1 #####
g1 <- gstat(NULL,"Comp1", Comp.1 ~ Macroarea + Macro_Period, individuals.spatial)
vm1 = variogram(g1, cutoff=1500)
plot(vm1, main="Variogram of PC1 in ILR space")

vm1.fit <- fit.variogram(vm1, vgm(0.2, "Exp", 500, 0.15))
plot(vm1, vm1.fit, pch = 3, main="Variogram of PC1 in ILR space")
vm1.fit

# plot with ggplot:
df.variogram1 <- data.frame(dist=vm1$dist, gamma=vm1$gamma)
preds <- variogramLine(vm1.fit, maxdist = max(df.variogram1$dist))

ggplot(df.variogram1, aes(x=dist, y=gamma)) +
  geom_point(size=3.5, colour='#3a9bdb') +
  geom_line(data=preds, colour='#3a9bdb', size=1.2) +
  ylim(0, 0.26) +
  theme(plot.title = element_text(hjust = 0.5), 
        panel.background = element_rect(fill = "#cae8ff",
                                        colour = "#cae8ff",
                                        size = 0.5, linetype = "solid")) +
  labs(y="Gamma", x="Geodetic distance") +
  ggtitle("Variogram of PC1 in ILR space")

# Create a gstat object setting an exponential variogram
# gstat(g.obj, id, formula, data, model, set,...)
g.tr1 <- gstat(formula = Comp.1 ~ Macroarea + Macro_Period, 
              data = individuals.spatial, model = vm1.fit)

# Let's check the coefficient:
coeff1 <- predict(g.tr1, all.comb, BLUE=T)
coeff1 <- data.frame(all.comb$Macroarea, all.comb$Macro_Period, coeff1$var1.pred)
coeff1[seq(6,80, by=12),]

##### Here we focus on the model for PC2 ######
g2 <- gstat(NULL,"Comp2", Comp.2 ~ Macroarea + Macro_Period, individuals.spatial)
vm2 = variogram(g2, cutoff=1500)
plot(vm2)

vm2.fit <- fit.variogram(vm2, vgm(0.07, "Exp", 1000, 0.05))
plot(vm2, vm2.fit, pch = 3, main="Variogram of PC2 in ILR space")
vm2.fit

# Create a gstat object setting an exponential variogram
# gstat(g.obj, id, formula, data, model, set,...)
g.tr2 <- gstat(formula = Comp.2 ~ Macroarea + Macro_Period, 
               data = individuals.spatial, model = vm2.fit)

# Let's check the coefficient:
coeff2 <- predict(g.tr2, all.comb, BLUE=T)
coeff2 <- data.frame(all.comb$Macroarea, all.comb$Macro_Period, coeff2$var1.pred)
coeff2[seq(6,80, by=12),]

##### Here we focus on the model for PC3 ######
g3 <- gstat(NULL,"Comp3", Comp.3 ~ Macroarea + Macro_Period, individuals.spatial)
vm3 = variogram(g3, cutoff=1500)
plot(vm3)

vm3.fit <- fit.variogram(vm3, vgm(0.01, "Exp", 1000, 0.015))
plot(vm3, vm3.fit, pch = 3, main="Variogram of PC3 in ILR space")
vm3.fit

# Create a gstat object setting an exponential variogram
# gstat(g.obj, id, formula, data, model, set,...)
g.tr3 <- gstat(formula = Comp.3 ~ Macroarea + Macro_Period, 
               data = individuals.spatial, model = vm3.fit)

# Let's check the coefficient:
coeff3 <- predict(g.tr3, all.comb, BLUE=T)
coeff3 <- data.frame(all.comb$Macroarea, all.comb$Macro_Period, coeff3$var1.pred)
coeff3[seq(6,80, by=12),]


##### Here we focus on the model for PC4 ######
g4 <- gstat(NULL,"Comp4", Comp.4 ~ Macroarea + Macro_Period, individuals.spatial)
vm4 = variogram(g4, cutoff=1500)
plot(vm4)

vm4.fit <- fit.variogram(vm4, vgm(0.0011, "Exp", 500, 0.008))
plot(vm4, vm4.fit, pch = 3, main="Variogram of PC4 in ILR space")
vm4.fit

# Create a gstat object setting an exponential variogram
# gstat(g.obj, id, formula, data, model, set,...)
g.tr4 <- gstat(formula = Comp.4 ~ Macroarea + Macro_Period, 
               data = individuals.spatial, model = vm4.fit)

# Let's check the coefficient:
coeff4 <- predict(g.tr4, all.comb, BLUE=T)
coeff4 <- data.frame(all.comb$Macroarea, all.comb$Macro_Period, coeff4$var1.pred)
coeff4[seq(6,80, by=12),]


#### Now we need to do the inverse CLR transformation #####

load = acomp(rbind(clrInv(pcx$loadings[,1]),
                   clrInv(pcx$loadings[,2]),
                   clrInv(pcx$loadings[,3]),
                   clrInv(pcx$loadings[,4])))

m.coeff <- matrix(0,84,5)

for (i in 1:84) {
  m.coeff[i,] <- coeff1[i,]$coeff1.var1.pred * acomp(load[1,]) + 
    coeff2[i,]$coeff2.var1.pred * acomp(load[2,]) + 
    coeff3[i,]$coeff3.var1.pred * acomp(load[3,]) +
    coeff4[i,]$coeff4.var1.pred * acomp(load[4,])
}

m.coeff <- as.data.frame(m.coeff)

m.coeff <- cbind(m.coeff, all.comb$Macroarea, all.comb$Macro_Period)
colnames(m.coeff) <- (c('WHG', 'EHG', 'Anatolia_N', 'CHG', 'North_Africa', 
                        'Macroarea', 'Macro_Period'))

m.coeff[which(m.coeff$Macroarea == 'Middle_East'), c(4,7,8)] 

offset <- length(unique(m.coeff$Macroarea))
m.coeff.ordered <- m.coeff
m.coeff.ordered[(offset+1):(2*offset),] <- m.coeff[(4*offset+1):(5*offset),]
m.coeff.ordered[(2*offset+1):(3*offset),] <- m.coeff[(3*offset+1):(4*offset),]
m.coeff.ordered[(3*offset+1):(4*offset),] <- m.coeff[(offset+1):(2*offset),]
m.coeff.ordered[(4*offset+1):(5*offset),] <- m.coeff[(6*offset+1):(7*offset),]
m.coeff.ordered[(5*offset+1):(6*offset),] <- m.coeff[(2*offset+1):(3*offset),]
m.coeff.ordered[(6*offset+1):(7*offset),] <- m.coeff[(5*offset+1):(6*offset),]
m.coeff.ordered <- m.coeff.ordered[,-1]
m.coeff.ordered

write.csv(m.coeff.ordered, 'data/comp_time_space.csv')

# plot of the composition of CHG in different Macroareas
quartz()
par(mfrow=c(3,4))
for (ma in unique(m.coeff$Macroarea)){
  barplot(height=m.coeff$CHG[m.coeff$Macroarea==ma], 
       ylab=paste('CHG composition', ma),
       names.arg=m.coeff$Macro_Period[m.coeff$Macroarea==ma],
       las=2, ylim=c(0.0, 0.3))
}