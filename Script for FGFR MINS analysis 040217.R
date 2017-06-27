#loadpackages
install.packages("dplyr")
install.packages("ggplot2")
install.packages("gplots")
install.packages("grid")
install.packages("plyr")
install.packages("reshape")
install.packages("reshape2")
install.packages("tidyr")
install.packages("xlsx")
install.packages("gridExtra")

#loadpackages
library("dplyr")
library("ggplot2")
library("gplots")
library("grid")
library("plyr")
library("reshape")
library("reshape2")
library("tidyr")
library("xlsx")
library("gridExtra")

setwd("~/Desktop/Z.FGFR.Rscript")
FR_all_clean <- read.csv("FR_all.csv")

#Convert values into Log-scale
FR_all_clean["logDAPI"]=log(FR_all_clean["CH1.Avg"])
FR_all_clean["logCDX2"]=log(FR_all_clean["CH3.Avg"])

#Correct CDX2 level based on Z using DAPI signal
lm(logDAPI~Z, data=FR_all_clean)
FR_all_clean["corrCDX2"]=FR_all_clean["logCDX2"]+FR_all_clean["Z"]*0.01764  ##Manually input slope ("Z" value) here

qplot(Z, corrCDX2, data=FR_all_clean)

#Plot cells with CDX2 and identify TE cells based on CDX2 level
FR_all_clean$Identity1 <- ifelse(FR_all_clean$corrCDX2>4, 'TE','ICM')

#Loading empirical bayse transformation as fuction
ebcor <- function(x, channel, group = NULL, logadd = 0.0001) {
  # get unique embryo IDs
  embryos <- unique(x$Embryo.Id)
  # if using integers (1 - 5) for channel, convert to standard names
  if (is.numeric(channel) == TRUE) {
    channels <- c('CH1.Avg', 'CH2.Avg', 'CH3.Avg', 'CH4.Avg', 'CH5.Avg')
    channel <- channels[channel]
  }
  # else, use whatever channel name given
  # data
  x$CH.Avg <- get(channel, x)
  # fitted regression coefficients and their standard errors
  coefs <- matrix(0, length(embryos), 2)
  for(i in 1:length(embryos)) {
    xi <- x[x$Embryo.Id==embryos[i],]
    coefs[i, 1:2] <- summary(lm(log(CH.Avg + logadd) ~ Z + (Identity1 == "TE"), 
                                data=xi))$coefficients[2, 1:2]
  }
  # if grouping variable is NULL create a dummy vector of 1s
  if (missing(group)) group <- rep(1, nrow(x))
  # group indicator of each embryo
  egrp <- tapply(group, x$Embryo.Id, function(y) {unique(y)[1]})
  # Emperical Bayes correction across the embryos in a group
  ebcoefs <- rep(0, length(embryos))
  for (i in unique(egrp)) {
    ebcoefs[egrp==i] <- mean(coefs[egrp==i, 1]) + 
      (1 - coefs[egrp==i, 2]^2/(coefs[egrp==i, 2]^2 + 
                                  var(coefs[egrp==i, 1])))*(coefs[egrp==i, 1] - 
                                                              mean(coefs[egrp==i, 1]))
  }
  # EB corrected log signal
  CH.ebLogCor <- rep(NA, nrow(x))
  for(i in 1:length(embryos)) {
    ii <- x$Embryo.Id==embryos[i]
    CH.ebLogCor[ii] <- log(logadd + x$CH.Avg[ii]) - ebcoefs[i]*x$Z[ii]
  }
  # return EB corrected log signal
  CH.ebLogCor
}


#Empirical baysian transformation of all channel
FR_all_clean["CH1.ebLogCor"] <- ebcor(FR_all_clean, channel=1)
FR_all_clean["CH2.ebLogCor"] <- ebcor(FR_all_clean, channel=2)
FR_all_clean["CH3.ebLogCor"] <- ebcor(FR_all_clean, channel=3)
FR_all_clean["CH4.ebLogCor"] <- ebcor(FR_all_clean, channel=4)
FR_all_clean["CH5.ebLogCor"] <- ebcor(FR_all_clean, channel=5)

#Convert corrected NAN-G6 values into lineage scale to assgin Identity
FR_all_clean$linCH4.ebLogCor <- exp(FR_all_clean$CH4.ebLogCor)
FR_all_clean$linCH5.ebLogCor <- exp(FR_all_clean$CH5.ebLogCor)

#Assign DP-EPI-PrE Identity based on NAN-G6 values
FR_all_clean$Identity3 <- ifelse(FR_all_clean$Identity1 == 'TE', 
                             'TE', 
                             ifelse((FR_all_clean$linCH4.ebLogCor < 450 &
                                       FR_all_clean$linCH5.ebLogCor < 450), 
                                    'DN', 
                                    ifelse((FR_all_clean$linCH4.ebLogCor > 450 &
                                              FR_all_clean$linCH5.ebLogCor < 450), 
                                           'EPI', 
                                           ifelse((FR_all_clean$linCH4.ebLogCor < 450 &
                                                     FR_all_clean$linCH5.ebLogCor > 450), 
                                                  'PrE', 'DP'))))


#factor identity
FR_all_clean$Identity3 <- factor(FR_all_clean$Identity3, levels=c('DN','DP','EPI','PrE','TE'))
FR_all_clean$Genotype <- as.factor(FR_all_clean$Genotype)

#Figure 3A-B
#Nanog-Gata6 expression in scatter plot
#logscale
ggplot(aes(x=CH5.ebLogCor, y=CH4.ebLogCor), data=subset(FR_all_clean, Identity1=="ICM"&Stage=="1"))+geom_point(aes(color=Identity3), size=0.5)+
  scale_color_manual(values=c(DN='darkgrey',EPI='red',DP='darkorchid2',PrE='cyan',TE='green'))+coord_fixed()+
  xlim(3,9)+ylim(3,9)+theme_bw()+facet_wrap(~Genotype)

ggplot(aes(x=CH5.ebLogCor, y=CH4.ebLogCor), data=subset(FR_all_clean, Identity1=="ICM"&Stage=="2"))+geom_point(aes(color=Identity3), size=0.5)+
  scale_color_manual(values=c(DN='darkgrey',EPI='red',DP='darkorchid2',PrE='cyan',TE='green'))+coord_fixed()+
  xlim(3,9)+ylim(3,9)+theme_bw()+facet_wrap(~Genotype)

ggplot(aes(x=CH5.ebLogCor, y=CH4.ebLogCor), data=subset(FR_all_clean, Identity1=="ICM"&Stage=="3"))+geom_point(aes(color=Identity3), size=0.5)+
  scale_color_manual(values=c(DN='darkgrey',EPI='red',DP='darkorchid2',PrE='cyan',TE='green'))+coord_fixed()+
  xlim(3,9)+ylim(3,9)+theme_bw()+facet_wrap(~Genotype)

#need this to use group_by function
detach(package:plyr)    
library(dplyr)

#Figure 3C-D
#Nanog or Gata6 levels in EPI or PrE
#Individual cell vs. AVG per embryo
#logscale
FR_all_clean <- FR_all_clean %>% group_by(Embryo.Id, Identity3) %>% mutate(AVGNG=mean(CH4.ebLogCor),AVGG6=mean(CH5.ebLogCor))
AVGFR <- select(FR_all_clean, Embryo.Id, Genotype, Stage, Identity3,AVGNG, AVGG6)
AVGFR <- as.data.frame(unique(AVGFR))
AVGFR$Genotype <- as.factor(AVGFR$Genotype)
AVGFR$Stage <- as.factor(AVGFR$Stage)

#Individual Cell
NG <- ggplot(aes(x=Genotype, y=CH4.ebLogCor), data=subset(FR_all_clean, Stage=="3"&Identity3=="EPI"))+
  facet_wrap(~Stage)+geom_boxplot(fill='red', outlier.shape=NA)+geom_jitter(size=0.1)+
  theme_bw()+scale_x_discrete(drop=F)

G6 <- ggplot(aes(x=Genotype, y=CH5.ebLogCor), data=subset(FR_all_clean, Stage=="3"&Identity3=="PrE"))+
  facet_wrap(~Stage)+geom_boxplot(fill='cyan',outlier.shape=NA)+geom_jitter(size=0.1)+
  theme_bw()+scale_x_discrete(drop=F)

grid.arrange(NG, G6, ncol=1)

#AVG per embryo
AVGNG <- ggplot(aes(x=Genotype, y=AVGNG), data=subset(AVGFR, Stage=="2"&Identity3=="EPI"))+
  facet_wrap(~Stage)+geom_boxplot(fill='red', outlier.shape = 1)+geom_jitter(size = 0.1)+
  theme_bw()+scale_x_discrete(drop=F)+ylim(6,9)
AVGG6 <- ggplot(aes(x=Genotype, y=AVGG6), data=subset(AVGFR, Stage=="2"&Identity3=="PrE"))+
  facet_wrap(~Stage)+geom_boxplot(fill='cyan', outlier.shape = 1)+geom_jitter(size = 0.1)+
  theme_bw()+scale_x_discrete(drop=F)+ylim(6,9)
grid.arrange(AVGNG, AVGG6, ncol=1)


####STATISTIACAL TEST WITH AVG NG and G6
#Statistical test - separate EPI vs.PrE
Stage1EPI <- subset(AVGFR, Stage=="1"&Identity3=="EPI")
Stage1PrE <- subset(AVGFR, Stage=="1"&Identity3=="PrE")
Stage2EPI <- subset(AVGFR, Stage=="2"&Identity3=="EPI")
Stage2PrE <- subset(AVGFR, Stage=="2"&Identity3=="PrE")
Stage3EPI <- subset(AVGFR, Stage=="3"&Identity3=="EPI")
Stage3PrE <- subset(AVGFR, Stage=="3"&Identity3=="PrE")

#Wilcox test compare WT vs. different genotypes
wilcox.test(AVGNG~Genotype, data=subset(Stage1EPI, Genotype=="1"|Genotype=="2"))
wilcox.test(AVGNG~Genotype, data=subset(Stage1EPI, Genotype=="1"|Genotype=="3"))
wilcox.test(AVGNG~Genotype, data=subset(Stage1EPI, Genotype=="1"|Genotype=="4"))
wilcox.test(AVGNG~Genotype, data=subset(Stage1EPI, Genotype=="1"|Genotype=="5"))
wilcox.test(AVGNG~Genotype, data=subset(Stage1EPI, Genotype=="1"|Genotype=="6"))
wilcox.test(AVGNG~Genotype, data=subset(Stage1EPI, Genotype=="1"|Genotype=="7"))
wilcox.test(AVGNG~Genotype, data=subset(Stage1EPI, Genotype=="1"|Genotype=="8"))
wilcox.test(AVGNG~Genotype, data=subset(Stage1EPI, Genotype=="1"|Genotype=="9"))

wilcox.test(AVGNG~Genotype, data=subset(Stage2EPI, Genotype=="1"|Genotype=="2"))
wilcox.test(AVGNG~Genotype, data=subset(Stage2EPI, Genotype=="1"|Genotype=="3"))
wilcox.test(AVGNG~Genotype, data=subset(Stage2EPI, Genotype=="1"|Genotype=="4"))
wilcox.test(AVGNG~Genotype, data=subset(Stage2EPI, Genotype=="1"|Genotype=="5"))
wilcox.test(AVGNG~Genotype, data=subset(Stage2EPI, Genotype=="1"|Genotype=="6"))
wilcox.test(AVGNG~Genotype, data=subset(Stage2EPI, Genotype=="1"|Genotype=="7"))
wilcox.test(AVGNG~Genotype, data=subset(Stage2EPI, Genotype=="1"|Genotype=="8"))
wilcox.test(AVGNG~Genotype, data=subset(Stage2EPI, Genotype=="1"|Genotype=="9"))

wilcox.test(AVGNG~Genotype, data=subset(Stage3EPI, Genotype=="1"|Genotype=="2"))
wilcox.test(AVGNG~Genotype, data=subset(Stage3EPI, Genotype=="1"|Genotype=="3"))
wilcox.test(AVGNG~Genotype, data=subset(Stage3EPI, Genotype=="1"|Genotype=="4"))
wilcox.test(AVGNG~Genotype, data=subset(Stage3EPI, Genotype=="1"|Genotype=="5"))
wilcox.test(AVGNG~Genotype, data=subset(Stage3EPI, Genotype=="1"|Genotype=="6"))
wilcox.test(AVGNG~Genotype, data=subset(Stage3EPI, Genotype=="1"|Genotype=="7"))
wilcox.test(AVGNG~Genotype, data=subset(Stage3EPI, Genotype=="1"|Genotype=="8"))
wilcox.test(AVGNG~Genotype, data=subset(Stage3EPI, Genotype=="1"|Genotype=="9"))

wilcox.test(AVGG6~Genotype, data=subset(Stage1PrE, Genotype=="1"|Genotype=="2"))
wilcox.test(AVGG6~Genotype, data=subset(Stage1PrE, Genotype=="1"|Genotype=="3"))
wilcox.test(AVGG6~Genotype, data=subset(Stage1PrE, Genotype=="1"|Genotype=="4"))
wilcox.test(AVGG6~Genotype, data=subset(Stage1PrE, Genotype=="1"|Genotype=="5"))
wilcox.test(AVGG6~Genotype, data=subset(Stage1PrE, Genotype=="1"|Genotype=="6"))
wilcox.test(AVGG6~Genotype, data=subset(Stage1PrE, Genotype=="1"|Genotype=="7"))
wilcox.test(AVGG6~Genotype, data=subset(Stage1PrE, Genotype=="1"|Genotype=="8"))
wilcox.test(AVGG6~Genotype, data=subset(Stage1PrE, Genotype=="1"|Genotype=="9"))

wilcox.test(AVGG6~Genotype, data=subset(Stage2PrE, Genotype=="1"|Genotype=="2"))
wilcox.test(AVGG6~Genotype, data=subset(Stage2PrE, Genotype=="1"|Genotype=="3"))
wilcox.test(AVGG6~Genotype, data=subset(Stage2PrE, Genotype=="1"|Genotype=="4"))
wilcox.test(AVGG6~Genotype, data=subset(Stage2PrE, Genotype=="1"|Genotype=="5"))
wilcox.test(AVGG6~Genotype, data=subset(Stage2PrE, Genotype=="1"|Genotype=="6"))
wilcox.test(AVGG6~Genotype, data=subset(Stage2PrE, Genotype=="1"|Genotype=="7"))
wilcox.test(AVGG6~Genotype, data=subset(Stage2PrE, Genotype=="1"|Genotype=="8"))
wilcox.test(AVGG6~Genotype, data=subset(Stage2PrE, Genotype=="1"|Genotype=="9"))

wilcox.test(AVGG6~Genotype, data=subset(Stage3PrE, Genotype=="1"|Genotype=="2"))
wilcox.test(AVGG6~Genotype, data=subset(Stage3PrE, Genotype=="1"|Genotype=="3"))
wilcox.test(AVGG6~Genotype, data=subset(Stage3PrE, Genotype=="1"|Genotype=="4"))
wilcox.test(AVGG6~Genotype, data=subset(Stage3PrE, Genotype=="1"|Genotype=="5"))
wilcox.test(AVGG6~Genotype, data=subset(Stage3PrE, Genotype=="1"|Genotype=="6"))
wilcox.test(AVGG6~Genotype, data=subset(Stage3PrE, Genotype=="1"|Genotype=="7"))
wilcox.test(AVGG6~Genotype, data=subset(Stage3PrE, Genotype=="1"|Genotype=="8"))
wilcox.test(AVGG6~Genotype, data=subset(Stage3PrE, Genotype=="1"|Genotype=="9"))

#When R2 level is reduced (het& KO), R1 affects NANOG levels in EPI
wilcox.test(AVGNG~Genotype, data=subset(Stage1EPI, Genotype=="4"|Genotype=="5"))
wilcox.test(AVGNG~Genotype, data=subset(Stage1EPI, Genotype=="7"|Genotype=="8"))
wilcox.test(AVGNG~Genotype, data=subset(Stage1EPI, Genotype=="8"|Genotype=="9"))
kruskal.test(AVGNG~Genotype, data=subset(Stage1EPI, Genotype=="4"|Genotype=="5"|Genotype=="6"))
kruskal.test(AVGNG~Genotype, data=subset(Stage1EPI, Genotype=="7"|Genotype=="8"|Genotype=="9"))

wilcox.test(AVGNG~Genotype, data=subset(Stage2EPI, Genotype=="4"|Genotype=="5"))
wilcox.test(AVGNG~Genotype, data=subset(Stage2EPI, Genotype=="7"|Genotype=="8"))
wilcox.test(AVGNG~Genotype, data=subset(Stage2EPI, Genotype=="8"|Genotype=="9"))
kruskal.test(AVGNG~Genotype, data=subset(Stage2EPI, Genotype=="4"|Genotype=="5"|Genotype=="6"))
kruskal.test(AVGNG~Genotype, data=subset(Stage2EPI, Genotype=="7"|Genotype=="8"|Genotype=="9"))

wilcox.test(AVGNG~Genotype, data=subset(Stage3EPI, Genotype=="4"|Genotype=="5"))
wilcox.test(AVGNG~Genotype, data=subset(Stage3EPI, Genotype=="7"|Genotype=="8"))
wilcox.test(AVGNG~Genotype, data=subset(Stage3EPI, Genotype=="8"|Genotype=="9"))
kruskal.test(AVGNG~Genotype, data=subset(Stage3EPI, Genotype=="4"|Genotype=="5"|Genotype=="6"))
kruskal.test(AVGNG~Genotype, data=subset(Stage3EPI, Genotype=="7"|Genotype=="8"|Genotype=="9"))
wilcox.test(AVGNG~Genotype, data=subset(Stage3EPI, Genotype=="6"|Genotype=="7"))


#When R2 level is reduced (het& KO), R1 affects GATA6 levels in PrE
wilcox.test(AVGG6~Genotype, data=subset(Stage1PrE, Genotype=="4"|Genotype=="5"))
wilcox.test(AVGG6~Genotype, data=subset(Stage1PrE, Genotype=="7"|Genotype=="8"))
wilcox.test(AVGG6~Genotype, data=subset(Stage1PrE, Genotype=="8"|Genotype=="9"))
kruskal.test(AVGG6~Genotype, data=subset(Stage1PrE, Genotype=="4"|Genotype=="5"|Genotype=="6"))
kruskal.test(AVGG6~Genotype, data=subset(Stage1PrE, Genotype=="7"|Genotype=="8"|Genotype=="9"))

wilcox.test(AVGG6~Genotype, data=subset(Stage2PrE, Genotype=="4"|Genotype=="5"))
wilcox.test(AVGG6~Genotype, data=subset(Stage2PrE, Genotype=="7"|Genotype=="8"))
wilcox.test(AVGG6~Genotype, data=subset(Stage2PrE, Genotype=="8"|Genotype=="9"))
kruskal.test(AVGG6~Genotype, data=subset(Stage2PrE, Genotype=="4"|Genotype=="5"|Genotype=="6"))
kruskal.test(AVGG6~Genotype, data=subset(Stage2PrE, Genotype=="7"|Genotype=="8"|Genotype=="9"))

wilcox.test(AVGG6~Genotype, data=subset(Stage3PrE, Genotype=="4"|Genotype=="5"))
wilcox.test(AVGG6~Genotype, data=subset(Stage3PrE, Genotype=="7"|Genotype=="8"))
wilcox.test(AVGG6~Genotype, data=subset(Stage3PrE, Genotype=="8"|Genotype=="9"))
kruskal.test(AVGG6~Genotype, data=subset(Stage3PrE, Genotype=="4"|Genotype=="5"|Genotype=="6"))
kruskal.test(AVGG6~Genotype, data=subset(Stage3PrE, Genotype=="7"|Genotype=="8"|Genotype=="9"))

sub <- dplyr::filter(FR_all_clean, Genotype == 9 & Stage == 3 & Identity3 != 'TE')
length(unique(sub$Embryo.Id))