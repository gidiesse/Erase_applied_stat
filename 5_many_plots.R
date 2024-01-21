# This is the file where we made most of the plots used to represent our data on the poster 

# loading the required packages
library(cowplot)      ## To work with ggpplot
library(ggplot2)      ## To make nice graphs
setwd("~/Desktop/ERASE_project2023_desgrp")

m.coeff <- read.csv('data/comp_time_space.csv')

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

## CHG over all periods of time - rainbow version - Central Europe
m.coeff.dat <- data.frame(CHG=m.coeff$CHG[m.coeff$Macroarea == 'Central_Europe'],
                          Time=factor(m.coeff$Macro_Period[m.coeff$Macroarea == 'Central_Europe'],
                                      levels=m.coeff$Macro_Period[m.coeff$Macroarea == 'Central_Europe']))

delta <- signif(m.coeff.dat$CHG[7] - m.coeff.dat$CHG[1],3) * 100

ggplot(m.coeff.dat, aes(x=Time, y=CHG, fill=Time)) + 
  scale_fill_manual(values = c('red', 'darkorange', 'yellow', 'green', 'dodgerblue', 
                               'mediumpurple2', 'hotpink')) +
  geom_bar(stat = "identity") + ggtitle('CHG Ancestry in Central Europe') + 
  theme(plot.title = element_text(hjust = 0.5), 
        panel.background = element_rect(fill = "#cae8ff",
                                        colour = "#cae8ff",
                                        size = 0.5, linetype = "solid")) +
  labs(fill = "Time Period") + geom_hline(yintercept=m.coeff.dat$CHG[7], 
                                          linetype="dashed", color = "hotpink") + 
  geom_hline(yintercept=m.coeff.dat$CHG[1], linetype="dashed", color = "red" ) +
  geom_segment(aes(x = 1, xend = 1, y = CHG[1], yend = CHG[7]), color = "black") + 
  geom_segment(aes(x = 0.9, xend = 1.1, y = CHG[1], yend = CHG[1]), color = "black") +
  geom_segment(aes(x = 0.9, xend = 1.1, y = CHG[7], yend = CHG[7]), color = "black") + 
  annotate("text", x=1.3, y=0.2, label= paste(delta, '%'))

## CHG over all periods of time for SICILY - chromatic scale 

# lightest green = #cff6c6
# light green = #69fe47
# medium green = #229b06
# dark green = #113808

m.coeff.dat <- data.frame(CHG=m.coeff$CHG[m.coeff$Macroarea == 'Sicily'],
                          Time=factor(m.coeff$Macro_Period[m.coeff$Macroarea == 'Sicily'],
                                      levels=m.coeff$Macro_Period[m.coeff$Macroarea == 'Sicily']))
m.coeff.dat <- data.frame(CHG=c(0.12, 0.17, 0.21, 0.26),
                          Time=factor(c(0,1,2,3)))

chg_sic <- ggplot(m.coeff.dat, aes(x=Time, y=CHG, fill=Time)) + 
  scale_fill_manual(values = c('#cff6c6', '#cff6c6', '#69fe47', '#69fe47', 
                               '#69fe47', '#69fe47', '#229b06'),
                    labels = c("5 - 10 %", "5 - 10 %", "10 - 15 % ", 
                               "10 - 15 % ", "10 - 15 % ", "10 - 15 % ", 
                               "15 - 20 %")) +
  geom_bar(stat = "identity") + ggtitle('CHG Ancestry in Sicily') + ylim(0.0, 0.3) +
  theme(plot.title = element_text(hjust = 0.5), 
        panel.background = element_rect(fill = "#cae8ff",
                                        colour = "#cae8ff", size = 0.5, linetype = "solid")) +
  geom_hline(yintercept=0.1, linetype="dashed", color = "black") +
  geom_hline(yintercept=0.15, linetype="dashed", color = "black") +
  geom_hline(yintercept=0.2, linetype="dashed", color = "black")

chg_sic <- ggplot(m.coeff.dat, aes(x=Time, y=CHG, fill=Time)) + 
  scale_fill_manual(values = c('#cff6c6','#69fe47', '#229b06', '#113808'),
                    labels = c("10 - 15 % ", "15 - 20 %", "20 - 25 %", "25 - 30 %"),
                    name = "CHG Composition") +
  geom_bar(stat = "identity") + ggtitle('CHG Ancestry in Sicily') + ylim(0.0, 0.3) +
  theme(plot.title = element_text(hjust = 0.5), 
        panel.background = element_rect(fill = "#cae8ff",
                                        colour = "#cae8ff", size = 0.5, linetype = "solid")) +
  geom_hline(yintercept=0.1, linetype="dashed", color = "black") +
  geom_hline(yintercept=0.15, linetype="dashed", color = "black") +
  geom_hline(yintercept=0.2, linetype="dashed", color = "black")

plot(chg_sic)

## CHG over all periods of time for EASTERN EUROPE - chromatic scale 
m.coeff.dat <- data.frame(CHG=m.coeff$CHG[m.coeff$Macroarea == 'Eastern_Europe'],
                          Time=factor(m.coeff$Macro_Period[m.coeff$Macroarea == 'Eastern_Europe'],
                                      levels=m.coeff$Macro_Period[m.coeff$Macroarea == 'Eastern_Europe']))

chg_asia <- ggplot(m.coeff.dat, aes(x=Time, y=CHG, fill=Time)) + 
  scale_fill_manual(values = c('#229b06', '#69fe47', '#113808', '#113808', 
                               '#113808', '#113808', '#113808'),
                    labels = c("20 - 25 %", "15 - 20 %", "25 - 30 % ", 
                               "25 - 30 % ", "25 - 30 % ", "25 - 30 % ", 
                               "25 - 30 %")) +
  geom_bar(stat = "identity") + ggtitle('CHG Ancestry in Eastern Europe') + ylim(0.0, 0.30) + 
  theme(plot.title = element_text(hjust = 0.5), 
        panel.background = element_rect(fill = "#cae8ff",
                                        colour = "#cae8ff", size = 0.5, linetype = "solid")) +
  geom_hline(yintercept=0.15, linetype="dashed", color = "black") +
  geom_hline(yintercept=0.2, linetype="dashed", color = "black") +
  geom_hline(yintercept=0.25, linetype="dashed", color = "black") 

plot(chg_asia)

plot_grid(chg_sic, chg_asia)


# Now we want two graphs to show how more isolated places have a lower change in 
# compared to "crossway" areas - we will compare Sardinia in the neolithic and in
# modern times with Balkans
our.anc <- c('WHG', 'EHG', 'Anatolia_N', 'CHG', 'North_Africa')

## Graph for Sardinia

idx.sard.neo <- which(m.coeff$Macroarea == 'Sardinia' &
                        m.coeff$Macro_Period == 'Neolithic')

m.coeff.dat1 <- data.frame(Sardinia=as.numeric(m.coeff[idx.sard.neo,2:6]),
                           Composition=factor(our.anc, levels=our.anc))

s_neo <- ggplot(m.coeff.dat1, aes(x=Composition, y=Sardinia, fill = Composition)) +
  ggtitle('Composition in Sardinia during Neolithic Period') +  
  theme(plot.title = element_text(hjust = 0.5), 
        panel.background = element_rect(fill = "#cae8ff",
                                        colour = "#cae8ff",
                                        size = 0.5, linetype = "solid")) +
  scale_fill_manual(values = c('red', 'darkorange', 'yellow', 'green', 'mediumblue')) +
  geom_bar(stat = "identity") + ylim(0.0, 0.5)

idx.sard.mod <- which(m.coeff$Macroarea == 'Sardinia' &
                        m.coeff$Macro_Period == 'ModernAge')


m.coeff.dat2 <- data.frame(Sardinia=as.numeric(m.coeff[idx.sard.mod,2:6]),
                           Composition=factor(our.anc, levels=our.anc))

s_mod <- ggplot(m.coeff.dat2, aes(x=Composition, y=Sardinia, fill = Composition)) +
  ggtitle('Composition in Sardinia during Modern Period') +  
  theme(plot.title = element_text(hjust = 0.5), 
        panel.background = element_rect(fill = "#cae8ff",
                                        colour = "#cae8ff",
                                        size = 0.5, linetype = "solid")) +
  scale_fill_manual(values = c('red', 'darkorange', 'yellow', 'green', 'mediumblue')) +
  geom_bar(stat = "identity") + ylim(0.0, 0.5)

plot_grid(s_neo, s_mod)

idx.sard.neo <- which(m.coeff$Macroarea == 'Sardinia' &
                        m.coeff$Macro_Period == 'Neolithic')

m.coeff.dat1 <- data.frame(Sardinia=as.numeric(m.coeff[idx.sard.neo,2:6]),
                           Composition=factor(our.anc, levels=our.anc))

s_neo <- ggplot(m.coeff.dat1, aes(x=Composition, y=Sardinia, fill = Composition)) +
  ggtitle('Composition in Sardinia during Neolithic Period') +  
  theme(plot.title = element_text(hjust = 0.5), 
        panel.background = element_rect(fill = "#cae8ff",
                                        colour = "#cae8ff",
                                        size = 0.5, linetype = "solid")) +
  scale_fill_manual(values = c('red', 'darkorange', 'yellow', 'green', 'mediumblue')) +
  geom_bar(stat = "identity") + ylim(0.0, 0.5)

idx.sard.mod <- which(m.coeff$Macroarea == 'Sardinia' &
                        m.coeff$Macro_Period == 'ModernAge')


m.coeff.dat2 <- data.frame(Sardinia=as.numeric(m.coeff[idx.sard.mod,2:6]),
                           Composition=factor(our.anc, levels=our.anc))

s_mod <- ggplot(m.coeff.dat2, aes(x=Composition, y=Sardinia, fill = Composition)) +
  ggtitle('Composition in Sardinia during Modern Period') +  
  theme(plot.title = element_text(hjust = 0.5), 
        panel.background = element_rect(fill = "#cae8ff",
                                        colour = "#cae8ff",
                                        size = 0.5, linetype = "solid")) +
  scale_fill_manual(values = c('red', 'darkorange', 'yellow', 'green', 'mediumblue')) +
  geom_bar(stat = "identity") + ylim(0.0, 0.5)

plot_grid(s_neo, s_mod)

## Graph for Balkans

idx.balk.neo <- which(m.coeff$Macroarea == 'Balkans' &
                        m.coeff$Macro_Period == 'Neolithic')

m.coeff.dat1 <- data.frame(Balkans=as.numeric(m.coeff[idx.balk.neo,2:6]),
                           Composition=factor(our.anc, levels=our.anc))

balk_neo <- ggplot(m.coeff.dat1, aes(x=Composition, y=Balkans, fill = Composition)) +
  ggtitle('Composition in the Balkans during Neolithic Period') +  
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c('red', 'darkorange', 'yellow', 'green', 'mediumblue')) +
  geom_bar(stat = "identity") + ylim(0.0, 0.5)

idx.balk.mod <- which(m.coeff$Macroarea == 'Balkans' &
                        m.coeff$Macro_Period == 'ModernAge')


m.coeff.dat2 <- data.frame(Balkans=as.numeric(m.coeff[idx.balk.mod,2:6]),
                           Composition=factor(our.anc, levels=our.anc))

balk_mod <- ggplot(m.coeff.dat2, aes(x=Composition, y=Balkans, fill = Composition)) +
  ggtitle('Composition in the Balkans during Modern Period') +  
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c('red', 'darkorange', 'yellow', 'green', 'mediumblue')) +
  geom_bar(stat = "identity") + ylim(0.0, 0.5)

plot_grid(balk_neo, balk_mod)

plot_grid(s_neo, s_mod, balk_neo, balk_mod)

## Graph for Western Europe

idx.we.neo <- which(m.coeff$Macroarea == 'Western_Europe' &
                      m.coeff$Macro_Period == 'Neolithic')

m.coeff.dat1 <- data.frame(WE=as.numeric(m.coeff[idx.we.neo,2:6]),
                           Composition=factor(our.anc, levels=our.anc))

we_neo <- ggplot(m.coeff.dat1, aes(x=Composition, y=WE, fill = Composition)) +
  ggtitle('Composition in Western Europe during Neolithic Period') +  
  theme(plot.title = element_text(hjust = 0.5), 
        panel.background = element_rect(fill = "#cae8ff",
                                        colour = "#cae8ff",
                                        size = 0.5, linetype = "solid")) +
  scale_fill_manual(values = c('red', 'darkorange', 'yellow', 'green', 'mediumblue')) +
  geom_bar(stat = "identity") + ylim(0.0, 0.5)

idx.we.mod <- which(m.coeff$Macroarea == 'Western_Europe' &
                      m.coeff$Macro_Period == 'ModernAge')


m.coeff.dat2 <- data.frame(WE=as.numeric(m.coeff[idx.we.mod,2:6]),
                           Composition=factor(our.anc, levels=our.anc))

we_mod <- ggplot(m.coeff.dat2, aes(x=Composition, y=WE, fill = Composition)) +
  ggtitle('Composition in Western Europe during Modern Period') +  
  theme(plot.title = element_text(hjust = 0.5), 
        panel.background = element_rect(fill = "#cae8ff",
                                        colour = "#cae8ff",
                                        size = 0.5, linetype = "solid")) +
  scale_fill_manual(values = c('red', 'darkorange', 'yellow', 'green', 'mediumblue')) +
  geom_bar(stat = "identity") + ylim(0.0, 0.5)

plot_grid(s_neo, s_mod, we_neo, we_mod)



################ MODERN PLOTS ###################

# 1) Northern Europe
idx.ne <- which(m.coeff$Macroarea == 'Northern_Europe' &
                  m.coeff$Macro_Period == 'ModernAge')

m.coeff.dat1 <- data.frame(NE=as.numeric(m.coeff[idx.ne,2:6]),
                           Composition=factor(our.anc, levels=our.anc))

ne_mod <- ggplot(m.coeff.dat1, aes(x=Composition, y=NE, fill = Composition)) +
  ggtitle('Composition in Northern Europe during Modern Age') +  
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c('red', 'darkorange', 'yellow', 'green', 'mediumblue')) +
  geom_bar(stat = "identity") + ylim(0.0, 0.4)

plot(ne_mod)

# 2) Central Europe
idx.ce <- which(m.coeff$Macroarea == 'Central_Europe' &
                  m.coeff$Macro_Period == 'ModernAge')

m.coeff.dat1 <- data.frame(CE=as.numeric(m.coeff[idx.ce,2:6]),
                           Composition=factor(our.anc, levels=our.anc))

ce_mod <- ggplot(m.coeff.dat1, aes(x=Composition, y=CE, fill = Composition)) +
  ggtitle('Composition in Central Europe during Modern Age') +  
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c('red', 'darkorange', 'yellow', 'green', 'mediumblue')) +
  geom_bar(stat = "identity") + ylim(0.0, 0.4)

plot(ce_mod)

# 3) Eastern Europe
idx.ee <- which(m.coeff$Macroarea == 'Eastern_Europe' &
                  m.coeff$Macro_Period == 'ModernAge')

m.coeff.dat1 <- data.frame(EE=as.numeric(m.coeff[idx.ee,2:6]),
                           Composition=factor(our.anc, levels=our.anc))

ee_mod <- ggplot(m.coeff.dat1, aes(x=Composition, y=EE, fill = Composition)) +
  ggtitle('Composition in Eastern Europe during Modern Age') +  
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c('red', 'darkorange', 'yellow', 'green', 'mediumblue')) +
  geom_bar(stat = "identity") + ylim(0.0, 0.4)

plot(ee_mod)

# 4) Western Europe
idx.we <- which(m.coeff$Macroarea == 'Western_Europe' &
                  m.coeff$Macro_Period == 'ModernAge')

m.coeff.dat1 <- data.frame(WE=as.numeric(m.coeff[idx.we,2:6]),
                           Composition=factor(our.anc, levels=our.anc))

we_mod <- ggplot(m.coeff.dat1, aes(x=Composition, y=WE, fill = Composition)) +
  ggtitle('Composition in Western Europe during Modern Age') +  
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c('red', 'darkorange', 'yellow', 'green', 'mediumblue')) +
  geom_bar(stat = "identity") + ylim(0.0, 0.4)

plot(we_mod)

# 5) Peninsular Italy
idx.pi <- which(m.coeff$Macroarea == 'Peninsular_Italy' &
                  m.coeff$Macro_Period == 'ModernAge')

m.coeff.dat1 <- data.frame(PI=as.numeric(m.coeff[idx.pi,2:6]),
                           Composition=factor(our.anc, levels=our.anc))

pi_mod <- ggplot(m.coeff.dat1, aes(x=Composition, y=PI, fill = Composition)) +
  ggtitle('Composition in Peninsular Italy during Modern Age') +  
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c('red', 'darkorange', 'yellow', 'green', 'mediumblue')) +
  geom_bar(stat = "identity") + ylim(0.0, 0.4)

plot(pi_mod)

# 6) Continental Italy
idx.ci <- which(m.coeff$Macroarea == 'Continental_Italy' &
                  m.coeff$Macro_Period == 'ModernAge')

m.coeff.dat1 <- data.frame(CI=as.numeric(m.coeff[idx.ci,2:6]),
                           Composition=factor(our.anc, levels=our.anc))

ci_mod <- ggplot(m.coeff.dat1, aes(x=Composition, y=CI, fill = Composition)) +
  ggtitle('Composition in Continental Italy during Modern Age') +  
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c('red', 'darkorange', 'yellow', 'green', 'mediumblue')) +
  geom_bar(stat = "identity") + ylim(0.0, 0.4)

plot(ci_mod)

# 7) Sardinia
idx.sard <- which(m.coeff$Macroarea == 'Sardinia' &
                    m.coeff$Macro_Period == 'ModernAge')

m.coeff.dat1 <- data.frame(Sard=as.numeric(m.coeff[idx.sard,2:6]),
                           Composition=factor(our.anc, levels=our.anc))

sard_mod <- ggplot(m.coeff.dat1, aes(x=Composition, y=Sard, fill = Composition)) +
  ggtitle('Composition in Sardinia during Modern Age') +  
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c('red', 'darkorange', 'yellow', 'green', 'mediumblue')) +
  geom_bar(stat = "identity") + ylim(0.0, 0.4)

plot(sard_mod)

# 8) Sicily
idx.sic <- which(m.coeff$Macroarea == 'Sicily' &
                   m.coeff$Macro_Period == 'ModernAge')

m.coeff.dat1 <- data.frame(Sic=as.numeric(m.coeff[idx.sic,2:6]),
                           Composition=factor(our.anc, levels=our.anc))

sic_mod <- ggplot(m.coeff.dat1, aes(x=Composition, y=Sic, fill = Composition)) +
  ggtitle('Composition in Sicily during Modern Age') +  
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c('red', 'darkorange', 'yellow', 'green', 'mediumblue')) +
  geom_bar(stat = "identity") + ylim(0.0, 0.4)

plot(sic_mod)

# 9) North Africa

idx.naf <- which(m.coeff$Macroarea == 'Northern_Africa' &
                   m.coeff$Macro_Period == 'ModernAge')

m.coeff.dat1 <- data.frame(N_africa=as.numeric(m.coeff[idx.naf,2:6]),
                           Composition=factor(our.anc, levels=our.anc))

naf_mod <- ggplot(m.coeff.dat1, aes(x=Composition, y=N_africa, fill = Composition)) +
  ggtitle('Composition in North Africa during Modern Age') +  
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c('red', 'darkorange', 'yellow', 'green', 'mediumblue')) +
  geom_bar(stat = "identity") + ylim(0.0, 0.4)

plot(naf_mod)

# 10) Balkans
idx.ba <- which(m.coeff$Macroarea == 'Balkans' &
                  m.coeff$Macro_Period == 'ModernAge')

m.coeff.dat1 <- data.frame(Balk=as.numeric(m.coeff[idx.ba,2:6]),
                           Composition=factor(our.anc, levels=our.anc))

ba_mod <- ggplot(m.coeff.dat1, aes(x=Composition, y=Balk, fill = Composition)) +
  ggtitle('Composition in the Balkans during Modern Age') +  
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c('red', 'darkorange', 'yellow', 'green', 'mediumblue')) +
  geom_bar(stat = "identity") + ylim(0.0, 0.4)

plot(ba_mod)

# 11) Middle East 
idx.me <- which(m.coeff$Macroarea == 'Middle_East' &
                  m.coeff$Macro_Period == 'ModernAge')

m.coeff.dat1 <- data.frame(ME=as.numeric(m.coeff[idx.me,2:6]),
                           Composition=factor(our.anc, levels=our.anc))

me_mod <- ggplot(m.coeff.dat1, aes(x=Composition, y=ME, fill = Composition)) +
  ggtitle('Composition in the Middle East during Modern Age') +  
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c('red', 'darkorange', 'yellow', 'green', 'mediumblue')) +
  geom_bar(stat = "identity") + ylim(0.0, 0.4)

plot(me_mod)

# 12) Asia
idx.as <- which(m.coeff$Macroarea == 'Asia' &
                  m.coeff$Macro_Period == 'ModernAge')

m.coeff.dat1 <- data.frame(Asia=as.numeric(m.coeff[idx.as,2:6]),
                           Composition=factor(our.anc, levels=our.anc))

me_as <- ggplot(m.coeff.dat1, aes(x=Composition, y=Asia, fill = Composition)) +
  ggtitle('Composition in Asia during Modern Age') +  
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c('red', 'darkorange', 'yellow', 'green', 'mediumblue')) +
  geom_bar(stat = "identity") + ylim(0.0, 0.4)

plot(me_as)

################ NEOLITHIC PLOTS ###################

# 1) Northern Europe
idx.ne <- which(m.coeff$Macroarea == 'Northern_Europe' &
                  m.coeff$Macro_Period == 'Neolithic')

m.coeff.dat1 <- data.frame(NE=as.numeric(m.coeff[idx.ne,2:6]),
                           Composition=factor(our.anc, levels=our.anc))

ne_mod <- ggplot(m.coeff.dat1, aes(x=Composition, y=NE, fill = Composition)) +
  ggtitle('Composition in Northern Europe during the Neolithic Period') +  
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c('red', 'darkorange', 'yellow', 'green', 'mediumblue')) +
  geom_bar(stat = "identity") + ylim(0.0, 0.5)

plot(ne_mod)

# 2) Central Europe
idx.ce <- which(m.coeff$Macroarea == 'Central_Europe' &
                  m.coeff$Macro_Period == 'Neolithic')

m.coeff.dat1 <- data.frame(CE=as.numeric(m.coeff[idx.ce,2:6]),
                           Composition=factor(our.anc, levels=our.anc))

ce_mod <- ggplot(m.coeff.dat1, aes(x=Composition, y=CE, fill = Composition)) +
  ggtitle('Composition in Central Europe during the Neolithic Period') +  
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c('red', 'darkorange', 'yellow', 'green', 'mediumblue')) +
  geom_bar(stat = "identity") + ylim(0.0, 0.5)

plot(ce_mod)

# 3) Eastern Europe
idx.ee <- which(m.coeff$Macroarea == 'Eastern_Europe' &
                  m.coeff$Macro_Period == 'Neolithic')

m.coeff.dat1 <- data.frame(EE=as.numeric(m.coeff[idx.ee,2:6]),
                           Composition=factor(our.anc, levels=our.anc))

ee_mod <- ggplot(m.coeff.dat1, aes(x=Composition, y=EE, fill = Composition)) +
  ggtitle('Composition in Eastern Europe during the Neolithic Period') +  
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c('red', 'darkorange', 'yellow', 'green', 'mediumblue')) +
  geom_bar(stat = "identity") + ylim(0.0, 0.5)

plot(ee_mod)

# 4) Western Europe
idx.we <- which(m.coeff$Macroarea == 'Western_Europe' &
                  m.coeff$Macro_Period == 'Neolithic')

m.coeff.dat1 <- data.frame(WE=as.numeric(m.coeff[idx.we,2:6]),
                           Composition=factor(our.anc, levels=our.anc))

we_mod <- ggplot(m.coeff.dat1, aes(x=Composition, y=WE, fill = Composition)) +
  ggtitle('Composition in Western Europe during the Neolithic Period') +  
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c('red', 'darkorange', 'yellow', 'green', 'mediumblue')) +
  geom_bar(stat = "identity") + ylim(0.0, 0.5)

plot(we_mod)

# 5) Peninsular Italy
idx.pi <- which(m.coeff$Macroarea == 'Peninsular_Italy' &
                  m.coeff$Macro_Period == 'Neolithic')

m.coeff.dat1 <- data.frame(PI=as.numeric(m.coeff[idx.pi,2:6]),
                           Composition=factor(our.anc, levels=our.anc))

pi_mod <- ggplot(m.coeff.dat1, aes(x=Composition, y=PI, fill = Composition)) +
  ggtitle('Composition in Peninsular Italy during the Neolithic Period') +  
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c('red', 'darkorange', 'yellow', 'green', 'mediumblue')) +
  geom_bar(stat = "identity") + ylim(0.0, 0.5)

plot(pi_mod)

# 6) Continental Italy
idx.ci <- which(m.coeff$Macroarea == 'Continental_Italy' &
                  m.coeff$Macro_Period == 'Neolithic')

m.coeff.dat1 <- data.frame(CI=as.numeric(m.coeff[idx.ci,2:6]),
                           Composition=factor(our.anc, levels=our.anc))

ci_mod <- ggplot(m.coeff.dat1, aes(x=Composition, y=CI, fill = Composition)) +
  ggtitle('Composition in Continental Italy during the Neolithic Period') +  
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c('red', 'darkorange', 'yellow', 'green', 'mediumblue')) +
  geom_bar(stat = "identity") + ylim(0.0, 0.5)

plot(ci_mod)

# 7) Sardinia
idx.sard <- which(m.coeff$Macroarea == 'Sardinia' &
                    m.coeff$Macro_Period == 'Neolithic')

m.coeff.dat1 <- data.frame(Sard=as.numeric(m.coeff[idx.sard,2:6]),
                           Composition=factor(our.anc, levels=our.anc))

sard_mod <- ggplot(m.coeff.dat1, aes(x=Composition, y=Sard, fill = Composition)) +
  ggtitle('Composition in Sardinia during the Neolithic Period') +  
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c('red', 'darkorange', 'yellow', 'green', 'mediumblue')) +
  geom_bar(stat = "identity") + ylim(0.0, 0.5)

plot(sard_mod)

# 8) Sicily
idx.sic <- which(m.coeff$Macroarea == 'Sicily' &
                   m.coeff$Macro_Period == 'Neolithic')

m.coeff.dat1 <- data.frame(Sic=as.numeric(m.coeff[idx.sic,2:6]),
                           Composition=factor(our.anc, levels=our.anc))

sic_mod <- ggplot(m.coeff.dat1, aes(x=Composition, y=Sic, fill = Composition)) +
  ggtitle('Composition in Sicily during the Neolithic Period') +  
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c('red', 'darkorange', 'yellow', 'green', 'mediumblue')) +
  geom_bar(stat = "identity") + ylim(0.0, 0.5)

plot(sic_mod)

# 9) North Africa

idx.naf <- which(m.coeff$Macroarea == 'Northern_Africa' &
                   m.coeff$Macro_Period == 'Neolithic')

m.coeff.dat1 <- data.frame(N_africa=as.numeric(m.coeff[idx.naf,2:6]),
                           Composition=factor(our.anc, levels=our.anc))

naf_mod <- ggplot(m.coeff.dat1, aes(x=Composition, y=N_africa, fill = Composition)) +
  ggtitle('Composition in North Africa during the Neolithic Period') +  
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c('red', 'darkorange', 'yellow', 'green', 'mediumblue')) +
  geom_bar(stat = "identity") + ylim(0.0, 0.5)

plot(naf_mod)

# 10) Balkans
idx.ba <- which(m.coeff$Macroarea == 'Balkans' &
                  m.coeff$Macro_Period == 'Neolithic')

m.coeff.dat1 <- data.frame(Balk=as.numeric(m.coeff[idx.ba,2:6]),
                           Composition=factor(our.anc, levels=our.anc))

ba_mod <- ggplot(m.coeff.dat1, aes(x=Composition, y=Balk, fill = Composition)) +
  ggtitle('Composition in the Balkans during the Neolithic Period') +  
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c('red', 'darkorange', 'yellow', 'green', 'mediumblue')) +
  geom_bar(stat = "identity") + ylim(0.0, 0.5)

plot(ba_mod)

# 11) Middle East 
idx.me <- which(m.coeff$Macroarea == 'Middle_East' &
                  m.coeff$Macro_Period == 'Neolithic')

m.coeff.dat1 <- data.frame(ME=as.numeric(m.coeff[idx.me,2:6]),
                           Composition=factor(our.anc, levels=our.anc))

me_mod <- ggplot(m.coeff.dat1, aes(x=Composition, y=ME, fill = Composition)) +
  ggtitle('Composition in the Middle East during the Neolithic Period') +  
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c('red', 'darkorange', 'yellow', 'green', 'mediumblue')) +
  geom_bar(stat = "identity") + ylim(0.0, 0.5)

plot(me_mod)

# 12) Asia
idx.as <- which(m.coeff$Macroarea == 'Asia' &
                  m.coeff$Macro_Period == 'Neolithic')

m.coeff.dat1 <- data.frame(Asia=as.numeric(m.coeff[idx.as,2:6]),
                           Composition=factor(our.anc, levels=our.anc))

me_as <- ggplot(m.coeff.dat1, aes(x=Composition, y=Asia, fill = Composition)) +
  ggtitle('Composition in Asia during the Neolithic Period') +  
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c('red', 'darkorange', 'yellow', 'green', 'mediumblue')) +
  geom_bar(stat = "identity") + ylim(0.0, 0.5)

plot(me_as)
