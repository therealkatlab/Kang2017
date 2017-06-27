#install packages
source("https://bioconductor.org/biocLite.R")
biocLite("pcaMethods")

#load packaged
library(dplyr)
library(plyr)
library(tidyr)
library(ggplot2)
library(gplots)
library(grid)
library(gridExtra)
library(reshape2)
library(RColorBrewer)
library(pcaMethods)

#Save combined csv file for all chips into one folder
##find folder manually
setwd("~/Desktop/Z.FGFR.Rscript")

read.data <- function(file){
  dat <- read.csv(file,header=TRUE)
  dat$fname <- file
  return(dat)
}

allChip <- do.call(rbind, lapply(list.files(pattern="csv$"),read.data))
rownames(allChip) <- allChip$X
allChip$X <- NULL
allChip$fname <- NULL

write.csv(allChip, "Fluidigm_chips.csv")

allChip <- read.csv("TableS6.csv", header = TRUE, stringsAsFactors = TRUE)
rownames(allChip) <- allChip$X
allChip$X <- NULL
allChip <- as.data.frame(t(allChip))
allChip$gene <- rownames(allChip)
gene <- read.csv("gene.csv", header = TRUE, stringsAsFactors = TRUE)
allChip <- dplyr::full_join(allChip, gene, by = "gene")
rownames(allChip) <- allChip$gene
write.csv(allChip, "TableS6_mod.csv")

#QC based on ct value (<8 or >28 failed)
allChip[(allChip)<8]<-28
allChip[(allChip)>=28]<-28

#calculate AG again (if one failed, replace with the other instead of average)
allChip$AG <- ifelse(allChip$Actb==28&allChip$Gapdh<28, allChip$Gapdh, ifelse(allChip$Gapdh==28&allChip$Actb<28, allChip$Actb, allChip$AG))

#Clean up the data set, assign EmbryoID & LitterID
allChip$CellID <- rownames(allChip)
allChip$EmbryoID <- as.factor(substr(allChip$CellID,1,2))
allChip$LitterID <- as.factor(substr(allChip$EmbryoID,1,1))

#Assign cross genotype for litter
allChip$Cross <- ifelse(allChip$LitterID=="B"|allChip$LitterID=="C"|allChip$LitterID=="E"|allChip$LitterID=="H"|
                          allChip$LitterID=="M"|allChip$LitterID=="O"|allChip$LitterID=="P"|allChip$LitterID=="Q"|allChip$LitterID=="S",'FR1','FR2')
allChip$Cross <- as.factor(allChip$Cross)

allChip$stage <- ifelse(allChip$LitterID=="A"|allChip$LitterID=="B"|
                          allChip$LitterID=="C"|allChip$LitterID=="D"|
                          allChip$LitterID=="E"|allChip$LitterID=="F"|allChip$LitterID=="G"|
                          allChip$LitterID=="H"|allChip$LitterID=="I"|
                          allChip$LitterID=="J"|allChip$LitterID=="K",'E3.5','E4.5')


#QC based on Actb, Gapdh and Sox2 levels [filter good quality cells and ICM cells]
allChip1 <- dplyr::filter(allChip, AG<20) #good quality
allChip2 <- dplyr::filter(allChip1, Sox2<20) # good quality ICM cells

#Normalize all Chip
##boxplot all cells per embryo
#Identify FGFR1 and FGFR2 Knockout embryos with normalized values
nallChip <- 28-allChip2[,1:49]
nallChip <- nallChip[,1:49]-nallChip[,49]
nallChip <- cbind(nallChip, allChip2[,50:54])

# 2 filters + normalized
allChip_clean <- nallChip

#Sort WT/FR1KO/FR2KO and recombine with genotype info
##WT embryos only include cells filtered after Fgfr1+QC+ICM
allChip_WT <- subset(allChip_clean, EmbryoID!="M2"&EmbryoID!="O3"&EmbryoID!="C6"&
                       EmbryoID!="C7"&EmbryoID!="E1"&EmbryoID!="M3"&EmbryoID!="O4"&
                       EmbryoID!="D4"&EmbryoID!="F5"&EmbryoID!="F1"&EmbryoID!="N4"&
                       EmbryoID!="N7"&EmbryoID!="T3"&EmbryoID!="N3"&EmbryoID!="N6"&EmbryoID!="U1"&EmbryoID!="V2"&
                       EmbryoID!="V3")
allChip_WT$Genotype <- 'WT'

##R1, R2 KO include all cells filtered after QC + ICM
allChip_clean_R1KO <- subset(allChip_clean, EmbryoID=="M2"|EmbryoID=="O3"|EmbryoID=="C6"|EmbryoID=="C7"|EmbryoID=="E1"|EmbryoID=="M3"|
                               EmbryoID=="O4")

allChip_clean_R1KO$Genotype <- 'R1KO'

allChip_clean_R2KO <- subset(allChip_clean, EmbryoID=="D4"|EmbryoID=="F5"|EmbryoID=="F1"|EmbryoID=="N4"|EmbryoID=="N7"|EmbryoID=="T3"|
                               EmbryoID=="N3"|EmbryoID=="N6"|EmbryoID=="U1"|EmbryoID=="V2"|EmbryoID=="V3")

allChip_clean_R2KO$Genotype <- 'R2KO'

allChip_clean <- rbind(allChip_WT, allChip_clean_R1KO, allChip_clean_R2KO)
allChip_backup <- allChip_clean ##backup the processed data as allChip_original (use this to reload)
#####################################################################################################

##################################Heatmap with all cells+all genes################################
allChip_clean$SI <- as.factor(paste(allChip_clean$stage, allChip_clean$Identity, sep=""))
allChip_clean$Identity2 <- ifelse(allChip_clean$Genotype=="R1KO", 'R1KO',ifelse(allChip_clean$Genotype=="R2KO", 'R2KO',allChip_clean$Identity))
allChip_clean$SI2 <- as.factor(paste(allChip_clean$stage, allChip_clean$Identity2, sep=""))

allChip_clean_heatmap <- allChip_clean[,1:46]## IT WAS CHANGED TO NOT TO INCLUDE AG, Actb, Gapdh
allChip_clean_heatmap <- as.matrix(allChip_clean_heatmap)
rownames(allChip_clean_heatmap) <- allChip_clean$SI2
allChip_clean_heatmap <- t(allChip_clean_heatmap)

hc.rows <- hclust(dist(allChip_clean_heatmap))
#plot(hc.rows)
hc.cols <- hclust(dist(t(allChip_clean_heatmap)))
#plot(hc.cols, cex=0.5)

#Heatmap with Z-score by row (by gene)
colors = c(seq(-2.5,0,length=10),seq(0,2.5,length=10))
colors=unique(colors)
my_palette <- colorRampPalette(c("#fff7bc","#fec44f","#d95f0e"))(n = 18) ## the n value needs to be unique color-1
#################################################################################################################
#Figure 5b
heatmap.2(allChip_clean_heatmap,trace='none',
          scale="row", 
          col=my_palette, breaks=colors)


########################Assign Identity with all 48 genes for WT + R1/R2 KOs############################
#Reload allChip_clean from allChip_backup
#Heatmap coloring may interfere PCA plots
#PCA performed here is only for lineage assignment purpose
allChip_clean <- allChip_backup
allChip_clean_WT <- subset(allChip_clean, Genotype=="WT")
allChip_clean_pca <- dplyr::select(allChip_clean_WT, -CellID, -EmbryoID, -LitterID, -Cross, -Genotype, -stage, -AG)
allChip_clean_pca <- as.matrix(allChip_clean_pca)
rownames(allChip_clean_pca) <- allChip_clean_WT$CellID

PCA = pca(data.matrix(allChip_clean_pca), method = "bpca", nPcs=3, center=TRUE, scale = "none")
PCA_score <- as.data.frame(PCA@scores)
allChip_clean_WT <- cbind(PCA_score, allChip_clean_WT)

###Lineage assignment with PCA scores+Nanog and GATA6 levles###
#PCA scores plot with NANOG-GATA6 levels to assign lineage identity
#Determine the cutoff lines for DP, EPI and PrE
G6<- ggplot(data=subset(allChip_clean_WT, Genotype=="WT"), aes(x=PC1, y=PC2, color=Gata6),target=row)+
  geom_point(size=2)+geom_density2d()+facet_wrap(~stage)
Nan <- ggplot(data=subset(allChip_clean_WT, Genotype=="WT"), aes(x=PC1, y=PC2, color=Nanog),target=row)+
  geom_point(size=2)+geom_density2d()+facet_wrap(~stage)
grid.arrange(G6, Nan, ncol=1)

#Separate E3.5 and E4.5
allChip_clean35 <- subset(allChip_clean_WT, stage=="E3.5")
allChip_clean45 <- subset(allChip_clean_WT, stage=="E4.5")

allChip_clean35$Identity <- ifelse(allChip_clean35$PC1<=0&abs(allChip_clean35$PC2)<(0.6)*(abs(allChip_clean35$PC1)), 'DP',
                                   ifelse(allChip_clean35$PC2>0, 'EPI','PrE'))
allChip_clean45$Identity <- ifelse(allChip_clean45$PC2>(-0.75), 'EPI','PrE')

qplot(PC1, PC2, data=allChip_clean35, color=Identity)
qplot(PC1, PC2, data=allChip_clean45, color=Identity)

allChip_clean_WT <- rbind(allChip_clean35, allChip_clean45)

#Separate R1KO and R2KO
allChip_clean_KO <- subset(allChip_clean, Genotype!="WT")
allChip_clean_pcaKO <- dplyr::select(allChip_clean_KO, -CellID, -EmbryoID, -LitterID, -Cross, -Genotype, -stage)
allChip_clean_pcaKO <- as.matrix(allChip_clean_pcaKO)
rownames(allChip_clean_pcaKO) <- allChip_clean_KO$CellID

PCAKO = pca(data.matrix(allChip_clean_pcaKO), method = "bpca", nPcs=3, center=TRUE, scale = "none")
PCAKO_score <- as.data.frame(PCAKO@scores)
allChip_clean_KO <- cbind(PCAKO_score, allChip_clean_KO)
PCAKO_score <- cbind(PCAKO_score, allChip_clean_KO)

allChip_cleanR1KO <- subset(allChip_clean_KO, Genotype=="R1KO")
allChip_cleanR2KO <- subset(allChip_clean_KO, Genotype=="R2KO")

allChip_cleanR1KO$Identity <- 'R1KO'
allChip_cleanR2KO$Identity <- 'R2KO'

allChip_clean_KO <- rbind(allChip_cleanR1KO, allChip_cleanR2KO)
allChip_clean <- rbind(allChip_clean_WT, allChip_clean_KO)

allChip_clean$SI <- as.factor(paste(allChip_clean$stage, allChip_clean$Identity, sep=""))
allChip_clean$Identity2 <- ifelse(allChip_clean$Genotype=="R1KO", 'R1KO',ifelse(allChip_clean$Genotype=="R2KO", 'R2KO',allChip_clean$Identity))
allChip_clean$SI2 <- as.factor(paste(allChip_clean$stage, allChip_clean$Identity2, sep=""))

######PCA analysis with 34 genes only with WT cells#######################
#1. Generate separate file with only 34 genes
allChip_clean_subset <- dplyr::select(allChip_clean,Pecam1, Esrrb, Fgf4, Nanog, Fgfr1, Etv5, Spry2, aPKCz,
                                      Etv4, Rex1, Sox2, Bmp4, Klf4, Pou5f1, Dab2, Stella, Spry4, Foxa2,Col4a1, Sox7,
                                      Zfp281, Bmpr2, Sox17, Pdgfra, Gata4, Lama1, Fgfr2, Xist,Stella,Grb2,
                                      Gata6, Dusp6, Actb, Gapdh, AG, EmbryoID, LitterID,CellID, 
                                      Genotype, stage, Identity, Identity2, SI, SI2) 

allChip_clean_subset_pca <- allChip_clean_subset[,c(-43,-42,-41,-40,-39,-38,-37,-36,-35)]
allChip_clean_subset_labels <- allChip_clean_subset[,c(43,42,41,40,39,38,37,36,35)]

#2. Generate file with only WT cells
WT_clean_subset <- subset(allChip_clean_subset, Genotype=="WT")

WT_clean_subset_pca <- WT_clean_subset[,c(-43,-42,-41,-40,-39,-38,-37,-36,-35)]
WT_clean_subset_labels <- WT_clean_subset[,c(43,42,41,40,39,38,37,36,35)]

PCA3 <- pca(data.matrix(WT_clean_subset_pca), method = "bpca", nPcs=3, center=TRUE, scale = "none")
PCA3_loadings <- as.data.frame(PCA3@loadings)
PCA3_score <- as.data.frame(PCA3@scores)

#3. Load prepared file with genes by their categories (for coloring)
gene <- read.csv("gene.csv")
PCA3_loadings$gene <- rownames(PCA3_loadings)
PCA3_loadings <- dplyr::left_join(PCA3_loadings, gene, by='gene')
PCA3_score <- cbind(PCA3_score, WT_clean_subset_pca, WT_clean_subset_labels)
###############################################################################################################

#Figure 5C
ggplot(data=PCA3_score, aes(x=PC1, y=PC2))+geom_point(size=2, aes(color=SI2))+theme_bw()+
  scale_color_manual(values=c('E3.5DP'='darkorchid1','E3.5EPI'='red','E3.5PrE'='cyan','E4.5DP'='darkorchid4',
                              'E4.5EPI'='firebrick','E4.5PrE'='cyan3','E3.5R1KO'='gold1','E4.5R1KO'='goldenrod3',
                              'E3.5R2KO'='gray70','E4.5R2KO'='gray33'))+coord_fixed()+ylim(-2.5,2.5)+xlim(-2.5,2.5)

#Figure 5D 
ggplot(PCA3_loadings, aes(x=PC1, y=PC2, label=gene, color=type), target=row)+
  geom_point(size=3, shape=15, aes(color=type))+geom_text(aes(color=type),size=2, hjust=0, nudge_x = 0.1)+
  coord_fixed(ratio=4/5)+labs(title="loadings WT only")+theme_bw()+xlim(-3.5,0.5)+ylim(-2,3)+
  scale_color_manual(values=c('EPI'='red','PrE'='cyan','FGF'='blue','others'='grey','H'='black', 'FGFR1'='goldenrod3', 'FGFR2'='gray33'))


######PCA analysis with 34 genes all cells#######################
#1. Generate separate file with only 34 genes
allChip_clean_subset <- dplyr::select(allChip_clean,Pecam1, Esrrb, Fgf4, Nanog, Fgfr1, Etv5, Spry2, aPKCz,
                                      Etv4, Rex1, Sox2, Bmp4, Klf4, Pou5f1, Dab2, Stella, Spry4, Col4a1, Sox7,
                                      Zfp281, Bmpr2, Sox17, Pdgfra, Gata4, Lama1, Fgfr2, Xist,Stella,Grb2,
                                      Foxa2, Gata6, Dusp6, Actb, Gapdh, AG, EmbryoID, LitterID,CellID, 
                                      Genotype, stage, Identity, Identity2, SI, SI2) 

allChip_clean_subset_pca <- allChip_clean_subset[,c(-43,-42,-42,-41,-40,-39,-38,-37,-36,-35)]
allChip_clean_subset_labels <- allChip_clean_subset[,c(43,42,41,40,39,38,37,36,35)]

PCA2 <- pca(data.matrix(allChip_clean_subset_pca), method = "bpca", nPcs=3, center=TRUE, scale = "none")

PCA2_loadings <- as.data.frame(PCA2@loadings)
PCA2_score <- as.data.frame(PCA2@scores)

## combine PCA with expression
PCA2_score <- cbind(PCA2_score, allChip_clean_subset_pca, allChip_clean_subset_labels)
PCA2_loadings$gene <- rownames(PCA2_loadings)
PCA2_loadings <- dplyr::left_join(PCA2_loadings, gene, by='gene')
################################################################################################################

#Figure 6B 
ggplot(data=PCA2_score, aes(x=PC1, y=PC2))+geom_point(size=2, aes(color=SI2))+theme_bw()+
  scale_color_manual(values=c('E3.5DP'='darkorchid1','E3.5EPI'='red','E3.5PrE'='cyan','E4.5DP'='darkorchid4',
                              'E4.5EPI'='firebrick','E4.5PrE'='cyan3','E3.5R1KO'='gold1','E4.5R1KO'='goldenrod3',
                              'E3.5R2KO'='gray70','E4.5R2KO'='gray33'))+coord_fixed()+ylim(-2.5,2.5)+xlim(-2.5,2.5)

#Figure 6C 
ggplot(PCA2_loadings, aes(x=PC1, y=PC2, label=gene), target=row)+
  geom_point(size=5, shape=15, aes(color=type))+geom_text(aes(color=type),size=5, hjust=0, nudge_x = 0.05)+
  coord_fixed(ratio=4/5)+labs(title="loadings all cells")+theme_bw()+xlim(-3.5,0.5)+ylim(-2,3)+
  scale_color_manual(values=c('EPI'='red','PrE'='cyan','FGF'='blue','others'='grey','H'='black', 'FGFR1'='goldenrod3', 'FGFR2'='gray33'))


######Analyze 34 genes by their groups
allChip_clean_subset_loadings <- dplyr::select(PCA2_score, -PC1,-PC2, -PC3,-EmbryoID, -LitterID)
allChip_clean_subset_loadings_ID <- dplyr::select(allChip_clean_subset_loadings, CellID, SI, stage, Genotype, Identity, SI2)
allChip_clean_subset_loadings_ID$CellID <- as.factor(as.character(allChip_clean_subset_loadings_ID$CellID))

rownames(allChip_clean_subset_loadings) <- allChip_clean_subset_loadings$CellID
allChip_clean_subset_loadings <- dplyr::select(allChip_clean_subset_loadings, -CellID, -SI, -stage, -Identity, -Genotype, -SI2, -Identity2)
allChip_clean_subset_loadings <- t(allChip_clean_subset_loadings)
allChip_clean_subset_loadings <- as.data.frame(allChip_clean_subset_loadings)

PCA2_loadings <- as.data.frame(PCA2@loadings)
PCA2_loadings <- cbind(PCA2_loadings, allChip_clean_subset_loadings)
PCA2_loadings$gene <- as.factor(rownames(PCA2_loadings))

PCA2_loadings <- melt(PCA2_loadings, id.vars=c('gene','PC1','PC2','PC3'))
PCA2_loadings$value <- as.numeric(as.character(PCA2_loadings$value))
colnames(PCA2_loadings)[colnames(PCA2_loadings)=='variable'] <- 'CellID'

PCA2_loadings <- dplyr::left_join(PCA2_loadings, gene, by='gene')
PCA2_loadings <- dplyr::left_join(PCA2_loadings, allChip_clean_subset_loadings_ID, by='CellID')

PCA2_loadings$type <- factor(PCA2_loadings$type, levels=c('EPI','PrE','FGF','others','H','FGFR1','FGFR2'))
PCA2_loadings$gene <- factor(PCA2_loadings$gene, levels = PCA2_loadings$gene[order(PCA2_loadings$type)])


AVGPCA2_loadings <- dplyr::select(PCA2_loadings, -CellID, -type, -SI, -SI2)
AVGPCA2_loadings$stage <- as.factor(AVGPCA2_loadings$stage)
AVGPCA2_loadings$Identity <- as.factor(AVGPCA2_loadings$Identity)
#AVGPCA2_loadings$GSI <- as.factor(paste(AVGPCA2_loadings$Genotype, AVGPCA2_loadings$SI, sep=""))
#AVGPCA2_loadings <- dplyr::select(AVGPCA2_loadings, -Genotype, -SI, -SI3)

library(dplyr)
detach(package:plyr)

AVGPCA2_loadings <- AVGPCA2_loadings %>% dplyr::group_by(stage, Genotype, Identity, gene) %>% summarise_each(funs(mean))
AVGPCA2_loadings <- dplyr::left_join(AVGPCA2_loadings, gene, by='gene')
AVGPCA2_loadings$type <- factor(AVGPCA2_loadings$type, levels=c('EPI','PrE','FGF','others','H', 'FGFR1','FGFR2'))
AVGPCA2_loadings$gene <- factor(AVGPCA2_loadings$gene, levels = AVGPCA2_loadings$gene[order(AVGPCA2_loadings$type)])
AVGPCA2_loadings$value2 <- 35+ AVGPCA2_loadings$value 
AVGPCA2_loadings$GI <- as.factor(paste(AVGPCA2_loadings$Genotype, AVGPCA2_loadings$Identity, sep=""))

#Figure 6D
ggplot(aes(x=type, y=value), data=subset(AVGPCA2_loadings, type!="FGFR1"&type!="FGFR2"))+geom_boxplot(aes(fill=type))+facet_grid(stage~GI)+
  scale_fill_manual(values=c('EPI'='red','PrE'='cyan','FGF'='blue','others'='grey','H'='black', 'FGFR1'='goldenrod3', 'FGFR2'='gray33'))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+coord_fixed(ratio=1/3) ##Figure 6C


plot1 <- ggplot(aes(x=type, y=value), data=subset(AVGPCA2_loadings, Genotype=="WT"&type!="FGFR1"&type!="FGFR2"))+geom_boxplot(aes(fill=type),size=0.3, outlier.size=0.3)+facet_grid(GI~stage)+
  scale_fill_manual(values=c('EPI'='red','PrE'='cyan','FGF'='blue','others'='grey','H'='black', 'FGFR1'='goldenrod3', 'FGFR2'='gray33'))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+coord_fixed(ratio=1/3)+ylim(-15,0.5)+theme_bw()+theme(legend.position="none") ##Figure 6C

plot2<- ggplot(aes(x=type, y=value), data=subset(AVGPCA2_loadings, Genotype=="R1KO"&type!="FGFR1"&type!="FGFR2"))+geom_boxplot(aes(fill=type),size=0.3, outlier.size=0.3)+facet_grid(GI~stage)+
  scale_fill_manual(values=c('EPI'='red','PrE'='cyan','FGF'='blue','others'='grey','H'='black', 'FGFR1'='goldenrod3', 'FGFR2'='gray33'))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+coord_fixed(ratio=1/3)+ylim(-15,0.5)+theme_bw()+theme(legend.position="none") ##Figure 6C

plot3<- ggplot(aes(x=type, y=value), data=subset(AVGPCA2_loadings, Genotype=="R2KO"&type!="FGFR1"&type!="FGFR2"))+geom_boxplot(aes(fill=type),size=0.3, outlier.size=0.3)+facet_grid(GI~stage)+
  scale_fill_manual(values=c('EPI'='red','PrE'='cyan','FGF'='blue','others'='grey','H'='black', 'FGFR1'='goldenrod3', 'FGFR2'='gray33'))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+coord_fixed(ratio=1/3)+ylim(-15,0.5)+theme_bw()+theme(legend.position="none") ##Figure 6C

grid.arrange(plot1, plot2, plot3, ncol=3)

PCA2_loadings$Genotype <- factor(PCA2_loadings$Genotype, levels=c('WT','R1KO','R2KO'))
PCA2_loadings$stage <- factor(PCA2_loadings$stage)

#Figure 6E
write.csv(PCA2_loadings, "PCA2_loadings.csv")
PCA2_loadings <- read.csv("PCA2_loadings.csv")

Spry2 <- ggplot(aes(x=SI2, y=value), data=subset(PCA2_loadings, gene=='Spry2'))+geom_boxplot(aes(fill=SI2),size=0.3,outlier.size=0.3)+
  scale_fill_manual(values=c('E3.5DP'='darkorchid1','E3.5EPI'='red','E3.5PrE'='cyan','E4.5DP'='darkorchid4',
                             'E4.5EPI'='firebrick','E4.5PrE'='cyan3','E3.5R1KO'='gold1','E4.5R1KO'='goldenrod3',
                             'E3.5R2KO'='gray70','E4.5R2KO'='gray33'))+theme_bw()+
  coord_fixed(ratio=1/3)+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
                               axis.ticks.x=element_blank(),strip.text.x = element_blank())+ theme(legend.position = "none")+labs(title="Spry2")+
  ylim(-15,0)+facet_wrap(~stage, scale="free_x")


Spry4 <- ggplot(aes(x=SI2, y=value), data=subset(PCA2_loadings, gene=='Spry4'))+geom_boxplot(aes(fill=SI2),size=0.3,outlier.size=0.3)+
  scale_fill_manual(values=c('E3.5DP'='darkorchid1','E3.5EPI'='red','E3.5PrE'='cyan','E4.5DP'='darkorchid4',
                             'E4.5EPI'='firebrick','E4.5PrE'='cyan3','E3.5R1KO'='gold1','E4.5R1KO'='goldenrod3',
                             'E3.5R2KO'='gray70','E4.5R2KO'='gray33'))+theme_bw()+
  coord_fixed(ratio=1/3)+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
                               axis.ticks.x=element_blank(),strip.text.x = element_blank())+ theme(legend.position = "none")+labs(title="Spry4")+
  ylim(-15,0)+facet_wrap(~stage, scale="free_x")



Etv4 <- ggplot(aes(x=SI2, y=value), data=subset(PCA2_loadings, gene=='Etv4'))+geom_boxplot(aes(fill=SI2),size=0.3,outlier.size=0.3)+
  scale_fill_manual(values=c('E3.5DP'='darkorchid1','E3.5EPI'='red','E3.5PrE'='cyan','E4.5DP'='darkorchid4',
                             'E4.5EPI'='firebrick','E4.5PrE'='cyan3','E3.5R1KO'='gold1','E4.5R1KO'='goldenrod3',
                             'E3.5R2KO'='gray70','E4.5R2KO'='gray33'))+theme_bw()+
  coord_fixed(ratio=1/3)+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
                               axis.ticks.x=element_blank(),strip.text.x = element_blank())+ theme(legend.position = "none")+labs(title="Etv4")+
  ylim(-15,0)+facet_wrap(~stage, scale="free_x")


Etv5 <- ggplot(aes(x=SI2, y=value), data=subset(PCA2_loadings, gene=='Etv5'))+geom_boxplot(aes(fill=SI2),size=0.3,outlier.size=0.3)+
  scale_fill_manual(values=c('E3.5DP'='darkorchid1','E3.5EPI'='red','E3.5PrE'='cyan','E4.5DP'='darkorchid4',
                             'E4.5EPI'='firebrick','E4.5PrE'='cyan3','E3.5R1KO'='gold1','E4.5R1KO'='goldenrod3',
                             'E3.5R2KO'='gray70','E4.5R2KO'='gray33'))+theme_bw()+
  coord_fixed(ratio=1/3)+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
                               axis.ticks.x=element_blank(),strip.text.x = element_blank())+ theme(legend.position = "none")+labs(title="Etv5")+
  ylim(-15,0)+facet_wrap(~stage, scale="free_x")


Dusp4 <- ggplot(aes(x=SI2, y=value), data=subset(PCA2_loadings, gene=='Dusp6'))+geom_boxplot(aes(fill=SI2),size=0.3,outlier.size=0.3)+
  scale_fill_manual(values=c('E3.5DP'='darkorchid1','E3.5EPI'='red','E3.5PrE'='cyan','E4.5DP'='darkorchid4',
                             'E4.5EPI'='firebrick','E4.5PrE'='cyan3','E3.5R1KO'='gold1','E4.5R1KO'='goldenrod3',
                             'E3.5R2KO'='gray70','E4.5R2KO'='gray33'))+theme_bw()+
  coord_fixed(ratio=1/3)+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
                               axis.ticks.x=element_blank(),strip.text.x = element_blank())+ theme(legend.position = "none")+labs(title="Dusp4")+
  ylim(-15,0)+facet_wrap(~stage, scale="free_x")


Gata6 <- ggplot(aes(x=SI2, y=value), data=subset(PCA2_loadings, gene=='Gata6'))+geom_boxplot(aes(fill=SI2),size=0.3,outlier.size=0.3)+
  scale_fill_manual(values=c('E3.5DP'='darkorchid1','E3.5EPI'='red','E3.5PrE'='cyan','E4.5DP'='darkorchid4',
                             'E4.5EPI'='firebrick','E4.5PrE'='cyan3','E3.5R1KO'='gold1','E4.5R1KO'='goldenrod3',
                             'E3.5R2KO'='gray70','E4.5R2KO'='gray33'))+theme_bw()+
  coord_fixed(ratio=1/3)+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
                               axis.ticks.x=element_blank(),strip.text.x = element_blank())+ theme(legend.position = "none")+labs(title="Gata6")+
  ylim(-15,0)+facet_wrap(~stage, scale="free_x")


Pou5f1 <- ggplot(aes(x=SI2, y=value), data=subset(PCA2_loadings, gene=='Pou5f1'))+geom_boxplot(aes(fill=SI2),size=0.3,outlier.size=0.3)+
  scale_fill_manual(values=c('E3.5DP'='darkorchid1','E3.5EPI'='red','E3.5PrE'='cyan','E4.5DP'='darkorchid4',
                             'E4.5EPI'='firebrick','E4.5PrE'='cyan3','E3.5R1KO'='gold1','E4.5R1KO'='goldenrod3',
                             'E3.5R2KO'='gray70','E4.5R2KO'='gray33'))+theme_bw()+
  coord_fixed(ratio=1/3)+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
                               axis.ticks.x=element_blank(),strip.text.x = element_blank())+ theme(legend.position = "none")+labs(title="Pou5f1")+
  ylim(-15,0)+facet_wrap(~stage, scale="free_x")

Dab2 <- ggplot(aes(x=SI2, y=value), data=subset(PCA2_loadings, gene=='Dab2'))+geom_boxplot(aes(fill=SI2),size=0.3,outlier.size=0.3)+
  scale_fill_manual(values=c('E3.5DP'='darkorchid1','E3.5EPI'='red','E3.5PrE'='cyan','E4.5DP'='darkorchid4',
                             'E4.5EPI'='firebrick','E4.5PrE'='cyan3','E3.5R1KO'='gold1','E4.5R1KO'='goldenrod3',
                             'E3.5R2KO'='gray70','E4.5R2KO'='gray33'))+theme_bw()+
  coord_fixed(ratio=1/3)+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
                               axis.ticks.x=element_blank(),strip.text.x = element_blank())+ theme(legend.position = "none")+labs(title="Dab2")+
  ylim(-15,0)+facet_wrap(~stage, scale="free_x")

Sox17 <- ggplot(aes(x=SI2, y=value), data=subset(PCA2_loadings, gene=='Sox17'))+geom_boxplot(aes(fill=SI2),size=0.3,outlier.size=0.3)+
  scale_fill_manual(values=c('E3.5DP'='darkorchid1','E3.5EPI'='red','E3.5PrE'='cyan','E4.5DP'='darkorchid4',
                             'E4.5EPI'='firebrick','E4.5PrE'='cyan3','E3.5R1KO'='gold1','E4.5R1KO'='goldenrod3',
                             'E3.5R2KO'='gray70','E4.5R2KO'='gray33'))+theme_bw()+
  coord_fixed(ratio=1/3)+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
                               axis.ticks.x=element_blank(),strip.text.x = element_blank())+ theme(legend.position = "none")+labs(title="Sox17")+
  ylim(-15,0)+facet_wrap(~stage, scale="free_x")

grid.arrange(Gata6,Dab2, Etv4, Pou5f1, Spry2, Etv5, Sox17, Spry4, Dusp4)

##Figure S6A
ggplot(aes(x=gene, y=value), data=subset(PCA2_loadings, stage=="E3.5"))+geom_boxplot(aes(fill=type))+facet_wrap(Genotype~Identity, ncol=3)+
  scale_fill_manual(values=c('EPI'='red','PrE'='cyan','FGF'='blue','others'='grey','H'='black', 'FGFR1'='goldenrod3', 'FGFR2'='gray33'))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
##Figure S6B
ggplot(aes(x=gene, y=value), data=subset(PCA2_loadings, stage=="E4.5"))+geom_boxplot(aes(fill=type))+facet_wrap(Genotype~Identity, ncol=2)+
  scale_fill_manual(values=c('EPI'='red','PrE'='cyan','FGF'='blue','others'='grey','H'='black', 'FGFR1'='goldenrod3', 'FGFR2'='gray33'))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


Gata6 <- ggplot(aes(x=SI2, y=value), data=subset(PCA2_loadings, gene=='Gata6'))+geom_boxplot(aes(fill=SI2),size=0.3,outlier.size=0.3)+
  scale_fill_manual(values=c('E3.5DP'='darkorchid1','E3.5EPI'='red','E3.5PrE'='cyan','E4.5DP'='darkorchid4',
                             'E4.5EPI'='firebrick','E4.5PrE'='cyan3','E3.5R1KO'='gold1','E4.5R1KO'='goldenrod3',
                             'E3.5R2KO'='gray70','E4.5R2KO'='gray33'))+theme_bw()+
  coord_fixed(ratio=1/3)+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
                               axis.ticks.x=element_blank(),strip.text.x = element_blank())+ theme(legend.position = "none")+labs(title="Gata6")+
  ylim(-15,0)+facet_wrap(~stage, scale="free_x")


Pdgfra <- ggplot(aes(x=SI2, y=value), data=subset(PCA2_loadings, gene=='Pdgfra'))+geom_boxplot(aes(fill=SI2),size=0.3,outlier.size=0.3)+
  scale_fill_manual(values=c('E3.5DP'='darkorchid1','E3.5EPI'='red','E3.5PrE'='cyan','E4.5DP'='darkorchid4',
                             'E4.5EPI'='firebrick','E4.5PrE'='cyan3','E3.5R1KO'='gold1','E4.5R1KO'='goldenrod3',
                             'E3.5R2KO'='gray70','E4.5R2KO'='gray33'))+theme_bw()+
  coord_fixed(ratio=1/3)+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
                               axis.ticks.x=element_blank(),strip.text.x = element_blank())+ theme(legend.position = "none")+labs(title="Pdgfra")+
  ylim(-15,0)+facet_wrap(~stage, scale="free_x")


Sox17 <- ggplot(aes(x=SI2, y=value), data=subset(PCA2_loadings, gene=='Sox17'))+geom_boxplot(aes(fill=SI2),size=0.3,outlier.size=0.3)+
  scale_fill_manual(values=c('E3.5DP'='darkorchid1','E3.5EPI'='red','E3.5PrE'='cyan','E4.5DP'='darkorchid4',
                             'E4.5EPI'='firebrick','E4.5PrE'='cyan3','E3.5R1KO'='gold1','E4.5R1KO'='goldenrod3',
                             'E3.5R2KO'='gray70','E4.5R2KO'='gray33'))+theme_bw()+
  coord_fixed(ratio=1/3)+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
                               axis.ticks.x=element_blank(),strip.text.x = element_blank())+ theme(legend.position = "none")+labs(title="Sox17")+
  ylim(-15,0)+facet_wrap(~stage, scale="free_x")


Gata4 <- ggplot(aes(x=SI2, y=value), data=subset(PCA2_loadings, gene=='Gata4'))+geom_boxplot(aes(fill=SI2),size=0.3,outlier.size=0.3)+
  scale_fill_manual(values=c('E3.5DP'='darkorchid1','E3.5EPI'='red','E3.5PrE'='cyan','E4.5DP'='darkorchid4',
                             'E4.5EPI'='firebrick','E4.5PrE'='cyan3','E3.5R1KO'='gold1','E4.5R1KO'='goldenrod3',
                             'E3.5R2KO'='gray70','E4.5R2KO'='gray33'))+theme_bw()+
  coord_fixed(ratio=1/3)+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
                               axis.ticks.x=element_blank(),strip.text.x = element_blank())+ theme(legend.position = "none")+labs(title="Gata4")+
  ylim(-15,0)+facet_wrap(~stage, scale="free_x")


Fgfr2 <- ggplot(aes(x=SI2, y=value), data=subset(PCA2_loadings, gene=='Sox7'))+geom_boxplot(aes(fill=SI2),size=0.3,outlier.size=0.3)+
  scale_fill_manual(values=c('E3.5DP'='darkorchid1','E3.5EPI'='red','E3.5PrE'='cyan','E4.5DP'='darkorchid4',
                             'E4.5EPI'='firebrick','E4.5PrE'='cyan3','E3.5R1KO'='gold1','E4.5R1KO'='goldenrod3',
                             'E3.5R2KO'='gray70','E4.5R2KO'='gray33'))+theme_bw()+
  coord_fixed(ratio=1/3)+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
                               axis.ticks.x=element_blank(),strip.text.x = element_blank())+ theme(legend.position = "none")+labs(title="Sox7")+
  ylim(-15,0)+facet_wrap(~stage, scale="free_x")

grid.arrange(Gata6,Pdgfra, Sox17,Gata4, Fgfr2)



Fgfr1 <- ggplot(aes(x=SI2, y=value), data=subset(PCA2_loadings, gene=='Fgfr1'))+geom_boxplot(aes(fill=SI2),size=0.3,outlier.size=0.3)+
  scale_fill_manual(values=c('E3.5DP'='darkorchid1','E3.5EPI'='red','E3.5PrE'='cyan','E4.5DP'='darkorchid4',
                             'E4.5EPI'='firebrick','E4.5PrE'='cyan3','E3.5R1KO'='gold1','E4.5R1KO'='goldenrod3',
                             'E3.5R2KO'='gray70','E4.5R2KO'='gray33'))+theme_bw()+ theme(legend.position = "none")+labs(title="Fgfr1")+
  ylim(-15,0)+facet_wrap(~Identity, scale="free_x")


Fgfr2 <- ggplot(aes(x=SI2, y=value), data=subset(PCA2_loadings, gene=='Fgfr2'))+geom_boxplot(aes(fill=SI2),size=0.3,outlier.size=0.3)+
  scale_fill_manual(values=c('E3.5DP'='darkorchid1','E3.5EPI'='red','E3.5PrE'='cyan','E4.5DP'='darkorchid4',
                             'E4.5EPI'='firebrick','E4.5PrE'='cyan3','E3.5R1KO'='gold1','E4.5R1KO'='goldenrod3',
                             'E3.5R2KO'='gray70','E4.5R2KO'='gray33'))+theme_bw()+
  coord_fixed(ratio=1/3)+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),
                               axis.ticks.x=element_blank(),strip.text.x = element_blank())+ theme(legend.position = "none")+labs(title="Fgfr2")+
  ylim(-15,0)+facet_wrap(~Identity, scale="free_x")

grid.arrange(Fgfr1,Fgfr2)

###FIGURE S5A
df <- dplyr::select(PCA2_loadings, -X, -PC1, -PC2, -PC3, -type)
df1 <- tidyr::spread(df, gene, value)

ggplot(aes(x = Gata6, y = Klf4), data = subset(df1, Genotype == "WT"))
+ geom_point(aes(fill = SI2), size=0.35) + coord_fixed() 
+ theme_bw() 
+ facet_wrap(~SI2)

ggplot(aes(x = Pou5f1, y = Klf4), data = subset(df1, Genotype == "WT")) + geom_point(aes(fill = SI2), size=0.35) + coord_fixed() + theme_bw() + facet_wrap(~SI2)

ggplot(aes(x = Sox2, y = Klf4), data = subset(df1, Genotype == "WT")) + geom_point(aes(fill = SI2), size=0.35) + coord_fixed() + theme_bw() + facet_wrap(~SI2)
