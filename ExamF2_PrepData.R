working_directory <-'/Users/Daisy/Desktop/MultivariateStats/multivar'
##Species Data####

#biomass data from Part D RDA
biomass <- read.delim(file = paste0(working_directory, '//data2_species.txt'), 
                      h=T, sep ="\t", row.names = 1)
biomass<-cast(biomass,siteID~species,sum)
abundance <-as.matrix(biomass)
#Species clusters from F.Part1####
clusters <-read.delim(file = paste0(working_directory, '//Examclusters.txt'), 
                      h=T, sep ="\t", row.names = 1)
clusters<- as.data.frame(t(clusters))
clusters<- clusters[-8]

##From "1.Exam F_Computing FDiv"
newdat <-read.delim(file = paste0(working_directory, '/Exam.SpeciesCoordinates.txt'), 
                    h=T, sep ="\t", row.names = 1)

#Temperature Data into 3 discrete categories####
env <-read.delim(file = paste0(working_directory, '//data2_env.txt'), 
                 h=T, sep ="\t", row.names = 1)
env.mean <- aggregate(env[, 3:5], list(env$siteID), mean)
## Load library and its dataset
require(dplyr)
boxplot(env.mean$Temperature)
T.quantile<-quantile(env.mean$Temperature,c(0,1/3,2/3,1)) #Quantiles
T.quantile[1]=T.quantile[1]-.00005
Temp <- env.mean %>% mutate(category=cut(Temperature, breaks=T.quantile,
                                         labels=c("low","moderate","high")))
boxplot(Temp$Temperature~Temp$category,col=3:5)
Temp <- data.frame(Temp[,-1], row.names = Temp[,1])
rm(env.mean,env)

#1.Low Temperature ######   
#subset of sites with low temperature
LowTemp <- Temp[Temp$category == "low", ]
#subset of species at sites with low Temperature (abundance = biomass.wide)
Species.Low<-as.data.frame(abundance[rownames(LowTemp),])
Species.Low["colsum",] <- colSums(Species.Low)
Species.Low["avg.den",] <- Species.Low["colsum",]/
  ((50*6)*(length(LowTemp$category)))
#50m2 per transect, 6 transects per site, #number of sites == individuals per m^2
AvgDen.Low<- Species.Low["avg.den",]

#Cluster 1 @Low Temp
Cluster1a <- which(clusters["x",]==1)
Cluster1.Low<-data.frame(AvgDen.Low[,Cluster1a])
rownames(Cluster1.Low) <- c("Cluster1.Low")

#Cluster 2 @Low Temp
Cluster2a <- which(clusters["x",]==2)
Cluster2.Low<-data.frame(AvgDen.Low[,Cluster2a])
rownames(Cluster2.Low) <- c("Cluster2.Low")

#2.Moderate Temperature ######   
#subset of sites with Moderate temperature
ModTemp <- Temp[Temp$category == "moderate", ]
#subset of species at sites with Mod Temperature (abundance = biomass.wide)
Species.Mod<-as.data.frame(abundance[rownames(ModTemp),])
Species.Mod["colsum",] <- colSums(Species.Mod)
Species.Mod["avg.den",] <- Species.Mod["colsum",]/
  ((50*6)*(length(ModTemp$category)))
AvgDen.Mod<- Species.Mod["avg.den",]

#Cluster 1 @Mod Temp
Cluster1b <- which(clusters["x",]==1)
Cluster1.Mod<-data.frame(AvgDen.Mod[,Cluster1b])
rownames(Cluster1.Mod) <- c("Cluster1.Mod")

#Cluster 2 @Mod Temp
Cluster2b <- which(clusters["x",]==2)
Cluster2.Mod<-data.frame(AvgDen.Mod[,Cluster2b])
rownames(Cluster2.Mod) <- c("Cluster2.Mod")

#3.High Temperature ######   
#subset of sites with High temperature

HighTemp <- Temp[Temp$category == "high", ]
#subset of species at sites with High Temperature (abundance = biomass.wide)
Species.High<-as.data.frame(abundance[rownames(HighTemp),])
Species.High["colsum",] <- colSums(Species.High)
Species.High["avg.den",] <- Species.High["colsum",]/
  ((50*6)*(length(HighTemp$category)))
AvgDen.High<- Species.High["avg.den",]

#Cluster 1 @High Temp
Cluster1c <- which(clusters["x",]==1)
Cluster1.High<-data.frame(AvgDen.High[,Cluster1c])
rownames(Cluster1.High) <- c("Cluster1.High")

#Cluster 2 @High Temp
Cluster2c <- which(clusters["x",]==2)
Cluster2.High<-data.frame(AvgDen.High[,Cluster2c])
rownames(Cluster2.High) <- c("Cluster2.High")

###Average Density
DenAvg.All<- rbind(AvgDen.Low,AvgDen.Mod,AvgDen.High)
rownames(DenAvg.All)<- c("AvgDen.Low","AvgDen.Mod","AvgDen.High")
write.table(DenAvg.All, file="DenAvg.All.txt",sep = "\t") 

##Cluster1####
Cluster1.All<- rbind(Cluster1.Low,Cluster1.Mod,Cluster1.High)
write.table(Cluster1.All, file="Cluster1.All.txt",sep = "\t") 

##Cluster2####
Cluster2.All<- rbind(Cluster2.Low,Cluster2.Mod,Cluster2.High)
write.table(Cluster2.All, file="Cluster2.All.txt",sep = "\t") 

