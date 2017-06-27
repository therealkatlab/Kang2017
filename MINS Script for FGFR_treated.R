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

FGFtreated_all <- read.csv("FGFtreated_all.csv")

###FGF4 1000ng/ml
FGFtreated <- subset(FGFtreated_all, Regime3=="FGF41000")

#Convert values into Log-scale
FGFtreated["logDAPI"]=log(FGFtreated["CH1.Avg"])
FGFtreated["logCDX2"]=log(FGFtreated["CH3.Avg"])
FGFtreated["logNG"]=log(FGFtreated["CH4.Avg"])
FGFtreated["logG6"]=log(FGFtreated["CH5.Avg"])
qplot(Z, logCDX2, data=FGFtreated)

#Correct CDX2 level based on Z using DAPI signal
lm(logDAPI~Z, data=subset(FGFtreated))
FGFtreated["corrCDX2"]=FGFtreated["logCDX2"]+FGFtreated["Z"]*0.01647

#Plot cells with CDX2 and identify TE cells based on CDX2 level
FGFtreated$Identity1 <- ifelse(FGFtreated$logCDX2>(-0.008)*FGFtreated$Z+3.5,'TE','ICM')

#Correct Nanog and Gata6 levels based on Z using CDX2 signal in TE cells
lm(logCDX2~Z, data=subset(FGFtreated, Identity1=="TE"))
FGFtreated["corrNG"]=FGFtreated["logNG"]+FGFtreated["Z"]*0.01497
FGFtreated["corrG6"]=FGFtreated["logG6"]+FGFtreated["Z"]*0.01497

FGFtreated$Identity2 <- ifelse(FGFtreated$Identity1=="TE",'TE',
                                  ifelse(FGFtreated$corrNG-FGFtreated$corrG6>0&FGFtreated$corrNG>4.5,'EPI',
                                         ifelse(FGFtreated$corrNG-(1.5)*(FGFtreated$corrG6)<(-3.75),'PrE',
                                                ifelse(FGFtreated$corrNG<4.5,'DN','DP'))))


FGFtreated$Identity2 <- factor(FGFtreated$Identity2, levels=c('DN','DP','EPI','PrE','TE'))
FGFtreated$Genotype <- as.factor(FGFtreated$Genotype)


ggplot(data=subset(FGFtreated, Identity1!="TE"), aes(Genotype, fill=Identity2))+
  geom_bar(position='fill')+
  scale_fill_manual(values=c(DN='darkgrey',EPI='red',DP='darkorchid2',PrE='cyan',TE='green')) # per genotype

ggplot(data=subset(FGFtreated, Identity1!="TE"), aes(Embryo.Id, fill=Identity2))+geom_bar(position='fill')+
  facet_wrap(~Genotype,ncol=3,scale="free_x")+
  scale_fill_manual(values=c(DN='darkgrey',EPI='red',DP='darkorchid2',PrE='cyan',TE='green'))# per embryo
write.csv(FGFtreated, "FGF4_1000.csv")

##FGF2 1000ng/ml
FGFtreated <- subset(FGFtreated_all, Regime3=="FGF21000")

#Convert values into Log-scale
FGFtreated["logDAPI"]=log(FGFtreated["CH1.Avg"])
FGFtreated["logCDX2"]=log(FGFtreated["CH3.Avg"])
FGFtreated["logNG"]=log(FGFtreated["CH4.Avg"])
FGFtreated["logG6"]=log(FGFtreated["CH5.Avg"])
qplot(Z, logCDX2, data=FGFtreated)

#Correct CDX2 level based on Z using DAPI signal
lm(logDAPI~Z, data=subset(FGFtreated))
FGFtreated["corrCDX2"]=FGFtreated["logCDX2"]+FGFtreated["Z"]*0.01429

#Plot cells with CDX2 and identify TE cells based on CDX2 level
qplot(Z, logCDX2, data=subset(FGFtreated))
FGFtreated$Identity1 <- ifelse(FGFtreated$logCDX2>(-0.01)*FGFtreated$Z+4,'TE','ICM')

#Correct Nanog and Gata6 levels based on Z using CDX2 signal in TE cells
lm(logCDX2~Z, data=subset(FGFtreated, Identity1=="TE"))
FGFtreated["corrNG"]=FGFtreated["logNG"]+FGFtreated["Z"]*0.01528
FGFtreated["corrG6"]=FGFtreated["logG6"]+FGFtreated["Z"]*0.01528


FGFtreated$Identity2 <- ifelse(FGFtreated$Identity1=="TE",'TE',
                                  ifelse((FGFtreated$corrNG-(5.5)*(FGFtreated$corrG6)<(-24.75)&FGFtreated$corrNG<5), 'PrE', 
                                         ifelse(FGFtreated$corrNG<5, 'DN',
                                                ifelse(FGFtreated$corrNG>FGFtreated$corrG6,'EPI','DP'))))



FGFtreated$Identity2 <- factor(FGFtreated$Identity2, levels=c('DN','DP','EPI','PrE','TE'))
FGFtreated$Genotype <- as.factor(FGFtreated$Genotype)


ggplot(data=subset(FGFtreated, Identity1!="TE"), aes(Genotype, fill=Identity2))+
  geom_bar(position='fill')+
  scale_fill_manual(values=c(DN='darkgrey',EPI='red',DP='darkorchid2',PrE='cyan',TE='green')) # per genotype

ggplot(data=subset(FGFtreated, Identity1!="TE"), aes(Embryo.Id, fill=Identity2))+geom_bar(position='fill')+
  facet_wrap(~Genotype,ncol=3,scale="free_x")+
  scale_fill_manual(values=c(DN='darkgrey',EPI='red',DP='darkorchid2',PrE='cyan',TE='green'))# per embryo
write.csv(FGFtreated, "FGF2_1000.csv")

###FGF4 2000ng/ml
FGFtreated <- subset(FGFtreated_all, Regime3=="FGF42000")

#Convert values into Log-scale
FGFtreated["logDAPI"]=log(FGFtreated["CH1.Avg"])
FGFtreated["logCDX2"]=log(FGFtreated["CH3.Avg"])
FGFtreated["logNG"]=log(FGFtreated["CH4.Avg"])
FGFtreated["logS17"]=log(FGFtreated["CH5.Avg"])
qplot(Z, logDAPI, data=FGFtreated)

#Correct CDX2 level based on Z using DAPI signal
lm(logDAPI~Z, data=subset(FGFtreated, logDAPI>0))
FGFtreated["corrCDX2"]=FGFtreated["logCDX2"]+FGFtreated["Z"]*0.0172

#Plot cells with CDX2 and identify TE cells based on CDX2 level
qplot(Z, logCDX2, data=subset(FGFtreated))
FGFtreated$Identity1 <- ifelse(FGFtreated$logCDX2>(-0.02)*FGFtreated$Z+5,'TE','ICM')

#Correct Nanog and Gata6 levels based on Z using CDX2 signal in TE cells
lm(logCDX2~Z, data=subset(FGFtreated, Identity1=="TE"))
FGFtreated["corrNG"]=FGFtreated["logNG"]+FGFtreated["Z"]*0.01947
FGFtreated["corrS17"]=FGFtreated["logS17"]+FGFtreated["Z"]*0.01947



FGFtreated$Identity2 <- ifelse(FGFtreated$Identity1=="TE",'TE',
                               ifelse(FGFtreated$corrNG>6.25&FGFtreated$corrS17>6,'DP',
                                      ifelse(FGFtreated$corrS17<5.75&FGFtreated$corrNG<5.3125,'DN',
                                             ifelse(FGFtreated$corrNG-3.75*(FGFtreated$corrS17)>(-16.25),'EPI','PrE'))))


FGFtreated$Identity2 <- factor(FGFtreated$Identity2, levels=c('DN','DP','EPI','PrE','TE'))
FGFtreated$Genotype <- as.factor(FGFtreated$Genotype)


ggplot(data=subset(FGFtreated, Identity1!="TE"), aes(Genotype, fill=Identity2))+
  geom_bar(position='fill')+
  scale_fill_manual(values=c(DN='darkgrey',EPI='red',DP='darkorchid2',PrE='cyan',TE='green')) # per genotype

ggplot(data=subset(FGFtreated, Identity1!="TE"), aes(Embryo.Id, fill=Identity2))+geom_bar(position='fill')+
  facet_wrap(~Genotype,ncol=3,scale="free_x")+
  scale_fill_manual(values=c(DN='darkgrey',EPI='red',DP='darkorchid2',PrE='cyan',TE='green'))# per embryo
write.csv(FGFtreated, "FGF4_2000.csv")