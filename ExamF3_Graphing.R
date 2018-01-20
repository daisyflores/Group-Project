working_directory <-'/Users/Daisy/Desktop/MultivariateStats/multivar'

#Average Density at sites with Low, Moderate, or High Temperature
Density <-read.delim(file = paste0(working_directory, '//DenAvg.All.txt'), 
                    h=T, sep ="\t", row.names = 1)
Density<-as.matrix(Density)
#Cluster1 
Cluster1 <- read.delim(file = paste0(working_directory, '//Cluster1.All.txt'), 
                       h=T, sep ="\t", row.names = 1)
Cluster1<-as.matrix(Cluster1)

#Cluster1 
Cluster2 <- read.delim(file = paste0(working_directory, '//Cluster2.All.txt'), 
                       h=T, sep ="\t", row.names = 1)
Cluster2<-as.matrix(Cluster2)

#Species Coordinates from F part 2
newdat <-read.delim(file = paste0(working_directory, '//Exam.SpeciesCoordinates.txt'), 
                      h=T, sep ="\t", row.names = 1)
newdat<- newdat[-8,]

# Global convex hull
vert.global<-data.frame(chull(newdat$A1,newdat$A2))
vert.global.coord<-newdat[c(vert.global[,1]),1:2]

# Low Temperature######                                                              
# Identity of species at low temperature 
sp_low<-which(Density["AvgDen.Low",]!=0) ##Species at sites with low T
set.low<-data.frame(newdat[sp_low,]) #Coordinates of those species
vert.low <- data.frame(chull(set.low$A1, set.low$A2)) # select what is going toConvex hull 
vert.low.coord<-set.low[c(vert.low[,1]),1:2] # Coordinates of the chull

# Species not feeding at low exposure
sp_0low<-which(Density["AvgDen.Low",]==0)
set.0low<-data.frame(newdat[sp_0low,])

# Species of Cluster 1 found at sites with low temperature
sp_clus1.low<-as.data.frame(which(Cluster1["Cluster1.Low",]!=0)) ## Cluster 1 
set.clus1.low<-data.frame(newdat[rownames(sp_clus1.low),1:2]) ## Cluster 1+ coordinates
vert.clus1.Low <- data.frame(chull(set.clus1.low$A1, set.clus1.low$A2))
vert.clus1.Low.coord<-set.clus1.low[c(vert.clus1.Low[,1]),1:2]

# Species of Cluster 2 found at sites with low temperature
sp_clus2.low<-as.data.frame(which(Cluster2["Cluster2.Low",]!=0))
set.clus2.low<-data.frame(newdat[rownames(sp_clus2.low),1:2]) #Cluster 2 + Coordinates
vert.clus2.Low <- data.frame(chull(set.clus2.low$A1, set.clus2.low$A2))
vert.clus2.Low.coord<-set.clus2.low[c(vert.clus2.Low[,1]),1:2]

# Relative frequency of bites at low exposure
nm_sp_low<-rownames(newdat[which(Density["AvgDen.Low",]>0),])
density_sp_low<-Density["AvgDen.Low",nm_sp_low]
rel.density_low<-density_sp_low/sum(density_sp_low)

## Plot
quartz() 
library(ggplot2)
plotLow<-ggplot(data=newdat,(aes(x=A1,y=A2)))+
  geom_polygon(data=vert.global,aes(x=vert.global.coord$A1,y=vert.global.coord$A2),
               fill=NA,color="lightgrey",lty=2)+
  geom_polygon(data=vert.low.coord,aes(x=A1,y=A2),
               fill="gray90",color="lightgrey",alpha=0.7,lwd=0.2)+
  geom_polygon(data=vert.clus1.Low.coord,aes(x=A1,y=A2),
               fill="#4D8963",color="#4D8963",alpha=0.3)+
  geom_polygon(data=vert.clus2.Low.coord,aes(x=A1,y=A2),
               fill="#41AAC4",color="#41AAC4",alpha=0.3)+
  geom_point(data=set.0low,aes(x=A1,y=A2),pch=4,cex=2)+
  geom_point(data=set.low,aes(x=A1,y=A2,size= rel.density_low),color="black",pch=16)+
  scale_size_continuous(name = "Intensity",limits=c(0,0.55),breaks = c(0.1,0.2,0.3,0.4,0.5), range = c(1,10))+
  scale_x_continuous(breaks = c(-3.0,-1.5,0,1.5,3.0))+
  scale_y_continuous(breaks = c(-3.0,-1.5,0,1.5))+
  guides(size=F)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
        title = element_text(vjust=1))+
  labs(title= "LOW",x = "PC 1", y = "PC 2")
plotLow #Plot Low####


#Moderate Temperature####
# Identity of species feeding at low exposure
sp_mod<-which(Density["AvgDen.Mod",]!=0)
set.mod<-data.frame(newdat[sp_mod,])
vert.mod <- data.frame(chull(set.mod$A1, set.mod$A2))
vert.mod.coord<-set.mod[c(vert.mod[,1]),1:2]

# Species not feeding at mod exposure
sp_0mod<-which(Density["AvgDen.Mod",]==0)
set.0mod<-data.frame(newdat[sp_0mod,])

# Species of Cluster 1 found at sites with mod temperature
sp_clus1.mod<-as.data.frame(which(Cluster1["Cluster1.Mod",]!=0)) ## Cluster 1 
set.clus1.mod<-data.frame(newdat[rownames(sp_clus1.mod),1:2]) ## Cluster 1+ coordinates
vert.clus1.Mod <- data.frame(chull(set.clus1.mod$A1, set.clus1.mod$A2))
vert.clus1.Mod.coord<-set.clus1.mod[c(vert.clus1.Mod[,1]),1:2]

# Species of Cluster 2 found at sites with mod temperature
sp_clus2.mod<-as.data.frame(which(Cluster2["Cluster2.Mod",]!=0))
set.clus2.mod<-data.frame(newdat[rownames(sp_clus2.mod),1:2]) #Cluster 2 + Coordinates
vert.clus2.Mod <- data.frame(chull(set.clus2.mod$A1, set.clus2.mod$A2))
vert.clus2.Mod.coord<-set.clus2.mod[c(vert.clus2.Mod[,1]),1:2]

# Relative frequency of bites at mod exposure
nm_sp_mod<-row.names(newdat)[which(Density["AvgDen.Mod",]>0)]
density_sp_mod<-Density["AvgDen.Mod",nm_sp_mod]
rel.density_mod<-density_sp_mod/sum(density_sp_mod) 

## Plot
quartz()
library(ggplot2)
plotmod<-ggplot(data=newdat,(aes(x=A1,y=A2)))+
  geom_polygon(data=vert.global,aes(x=vert.global.coord$A1,y=vert.global.coord$A2),
               fill=NA,color="lightgrey",lty=2)+
  geom_polygon(data=vert.mod.coord,aes(x=A1,y=A2),
               fill="gray90",color="lightgrey",alpha=0.7,lwd=0.2)+
  geom_polygon(data=vert.clus1.Mod.coord,aes(x=A1,y=A2),
               fill="#4D8963",color="#4D8963",alpha=0.3)+
  geom_polygon(data=vert.clus2.Mod.coord,aes(x=A1,y=A2),
               fill="#41AAC4",color="#41AAC4",alpha=0.3)+
  geom_point(data=set.0mod,aes(x=A1,y=A2),pch=4,cex=2)+
  geom_point(data=set.mod,aes(x=A1,y=A2,size=rel.density_mod),color="black",pch=16)+
  scale_size_continuous(name = "Intensity",limits=c(0,0.55),breaks = c(0.1,0.2,0.3,0.4,0.5), range = c(1,10))+
  scale_x_continuous(breaks = c(-3.0,-1.5,0,1.5,3.0))+
  scale_y_continuous(breaks = c(-3.0,-1.5,0,1.5))+
  guides(size=F)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
        title = element_text(vjust=1))+
  labs(title= "MOD",x = "PC 1", y = "PC 2")
plotmod #Plot mod####

#High Temperature####
sp_high<-which(Density["AvgDen.High",]!=0)
set.high<-data.frame(newdat[sp_high,])
vert.high <- data.frame(chull(set.high$A1, set.high$A2))
vert.high.coord<-set.high[c(vert.high[,1]),1:2]

# Species not feeding at high exposure
sp_0high<-which(Density["AvgDen.High",]==0)
set.0high<-data.frame(newdat[sp_0high,])

# Species of Cluster 1 found at sites with high temperature
sp_clus1.high<-as.data.frame(which(Cluster1["Cluster1.High",]!=0)) ## Cluster 1 
set.clus1.high<-data.frame(newdat[rownames(sp_clus1.high),1:2]) ## Cluster 1+ coordinates
vert.clus1.High <- data.frame(chull(set.clus1.high$A1, set.clus1.high$A2))
vert.clus1.High.coord<-set.clus1.high[c(vert.clus1.High[,1]),1:2]

# Species of Cluster 2 found at sites with high temperature
sp_clus2.high<-as.data.frame(which(Cluster2["Cluster2.High",]!=0))
set.clus2.high<-data.frame(newdat[rownames(sp_clus2.high),1:2]) #Cluster 2 + Coordinates
vert.clus2.High <- data.frame(chull(set.clus2.high$A1, set.clus2.high$A2))
vert.clus2.High.coord<-set.clus2.high[c(vert.clus2.High[,1]),1:2]

# Relative frequency of bites at high exposure
nm_sp_high<-row.names(newdat)[which(Density["AvgDen.High",]>0)]
density_sp_high<-Density["AvgDen.High",nm_sp_high]
rel.density_high<-density_sp_high/sum(density_sp_high) 

## Plot
quartz()
library(ggplot2)
plothigh<-ggplot(data=newdat,(aes(x=A1,y=A2)))+
  geom_polygon(data=vert.global,aes(x=vert.global.coord$A1,y=vert.global.coord$A2),
               fill=NA,color="lightgrey",lty=2)+
  geom_polygon(data=vert.high.coord,aes(x=A1,y=A2),
               fill="gray90",color="lightgrey",alpha=0.7,lwd=0.2)+
  geom_polygon(data=vert.clus1.High.coord,aes(x=A1,y=A2),
               fill="#4D8963",color="#4D8963",alpha=0.3)+
  geom_polygon(data=vert.clus2.High.coord,aes(x=A1,y=A2),
               fill="#41AAC4",color="#41AAC4",alpha=0.3)+
  geom_point(data=set.0high,aes(x=A1,y=A2),pch=4,cex=2)+
  geom_point(data=set.high,aes(x=A1,y=A2,size=rel.density_high),color="black",pch=16)+
  scale_size_continuous(name = "Intensity",limits=c(0,0.55),breaks = c(0.1,0.2,0.3,0.4,0.5), range = c(1,10))+
  scale_x_continuous(breaks = c(-3.0,-1.5,0,1.5,3.0))+
  scale_y_continuous(breaks = c(-3.0,-1.5,0,1.5))+
  guides(size=F)+
  theme_bw()+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
        title = element_text(vjust=1))+
  labs(title= "HIGH",x = "PC 1", y = "PC 2")
plothigh #plot high####
dev.off()
quartz()
par(mfrow= c(1,3))

