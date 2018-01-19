rm(list = ls())
working_directory <-'/Users/Daisy/Desktop/MultivariateStats/multivar'

# 1.1. ORDINATION OF SPECIES BASED ON MORPHOLOGICAL TRAITS 

#1. Compute FD using Multidim FD####
source("quality_funct_space.R")
source("plot_funct_space.R")
source("multidimFD.R")
library(ggplot2)
library(grid)
library(vegan)
library(ade4)
library(ggthemes)
library(NbClust)

# Import traits matrix
traits<- read.table(paste0(working_directory, '//urchingtraits.txt'), 
                    h=T, row.names = 1)

# A. RUN THE ORDINATION
traits<-scale(traits)
d<-dist(traits,method = "euclidean")
pcoa<-dudi.pco(d,scannf = F,full = T)

traitspca<-rda(traits, scale=TRUE)
traitspca
summary(traitspca, scaling=1)
plot(traitspca)
cleanplot.pca(traitspca)

# B. EXTRACT EIGEN VALUES 
eig<-pcoa$eig 
eig/sum(eig)# first captures 50% and than 20%
cumsum(eig/sum(eig))# comulative percentages
# How much variability is captured by the axes?
# How many axes should be retained?

# C. PLOTTING 
# EXTRACT SCORES FOR SPECIES AND VECTORS (in this case vectors = traits)
scrs<- as.data.frame(scores(pcoa, display = "sites"))
vf <- envfit(pcoa, traits, perm = 999)
spp.scrs <- as.data.frame(scores(vf, display = "vectors"))

write.table(scrs, file="Exam.SpeciesCoordinates.txt", sep="\t") 

# PLOTTING VECTORS FIRST
# CREATE FULL LENGTH NAMES OF TRAITS
spp.scrs <- cbind(spp.scrs, 
                  Traits=c("Maximum size","Erosion rate","Sedimentation rate",
                           "Gut to jaw ratio", "Jaw size", "Div gut micro"))
VectorPlot<-ggplot(data = scrs, aes(A1, A2)) + 
  geom_segment(data=spp.scrs,aes(x=0,xend=A1,y=0,yend=A2),
               arrow = arrow(length = unit(0, "cm")),colour="darkgrey") + 
  geom_text(data=spp.scrs,aes(x=A1+0.06,y=A2,label=Traits),size=6)+
  xlim(-1,1.3)+
  ylim(-1,1.3)+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text = element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_blank())
VectorPlot

# 1.2. WARD'S CLUSTERING OF SPECIES AND PLOTTING SPECIES ORDINATION WITH CLUSTERING INFORMATION
# OVERLAYED
par(mfrow=c(1,1)) 
clust <- hclust(d, method="ward.D2") #clustering analysis

# Determining the best number of clusters 
# IN THIS OCCASION WE USED THE PACKAGE NbClusT TO SELECT THE BEST NUMBER OF CLUSTERS
# THIS PACKAGE HANDLES THE AMBIGUITY OF THE DIFFERENT INDEXES AND COMPUTES
# MANY INDICES, IT THEN DECIDES ON THE BEST NUMBER OF CLUSTERS BASED ON WHAT 
# THE MAJORITY OF INDICES INDICATED 
library(NbClust)
nb <- NbClust(traits,distance="euclidean", min.nc=2, max.nc=5, 
              method = "ward.D2", index = "all",alphaBeale = 0.01)
par(mfrow=c(1,1))
hist(nb$Best.nc[1,], breaks = max(na.omit(nb$Best.nc[1,])), 
     xlab="Number of clusters",col = "lightgrey", main='')
# Beales method tells us that the most sound decision is to cut the tree in 2 
clusters <- cutree(clust, k=2)
clusters<-as.factor(clusters) 
write.table(clusters, file="Examclusters.txt", sep="\t") 


# extracting the species names for plotting
speciesnames<-rownames(scrs)

# FINAL ORDINATION PLOT ###############################################
OrdiPLot<-ggplot(data = scrs, aes(A1, A2,color=clusters))+
  geom_point(aes(x=A1, y=A2,shape=clusters),size=3.5)+
  scale_shape_manual(values=c(15,17),guide=F)+
  geom_text(aes(x=A1, y=A2+0.15,label=speciesnames),size=3)+
  scale_color_manual(values=c("black","orange"),guide=FALSE)+
  labs(x="PC 1 (49.9 %)",y="PC 2 (20.3 %)")+
  theme(legend.position="none",
        axis.text.x = element_text(size=13),
        axis.text.y = element_text(size=13),
        axis.title.x=element_text(size=15,vjust=-0.3),
        axis.title.y=element_text(size=15,vjust=1),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
OrdiPLot
# x axis limits
OrdiPLot + xlim(-3, 2)
# y axis limits
OrdiPLot + ylim(-2.5, 2)
# USE GGPLOT THEME ASSISTANT TO IMPROVE THE PLOT #######################

# THE END ##############################################################

# SO ~70% OF THE SPECIES CAN BE DEVIDID INTO 2 MAJOR CLUSTERS DRIVEN BY
# THE TRAITS/FUNCTION


