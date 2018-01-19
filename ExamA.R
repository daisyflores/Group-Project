working_directory <-'/Users/Daisy/Desktop/MultivariateStats/multivar'

#1.Read Data####
env <-read.delim(file = paste0(working_directory, '//data2_env.txt'), 
                 h=T, sep ="\t", row.names = 1)
species <-(read.delim(file = paste0(working_directory, '//data2_species.txt'), 
                      h=T, sep ="\t"))

#2.Calculate density/site ####
require(reshape) # (load necessary package)
species.wide<-cast(species,siteID+transect~species,sum)
species.melted <- melt(species.wide)
row.names(species.melted) <- 1:length(species.melted[,1])
names(species.melted)[3] <- c("abundance")
#disregards zeros
species.density<- aggregate(abundance~ siteID,
                            FUN=sum,data=species.melted)
species.density$density <- species.density$abundance/(50*6) #50m2 * 6 transects

##Prepare to combine into one table ####
latitude<- aggregate(Latitude~ siteID,
                            FUN=mean,env)


all.equal(as.character(latitude$siteID),
          as.character(species.density$siteID))

#Combine Data
distribution <- data.frame(species.density,latitude$Latitude)
row.names(distribution) <- 1:length(distribution[,1])
names(distribution)[ncol(distribution)] <- "latitude"
distribution.ordered <- distribution[order(distribution$latitude),]

write.table(distribution, file = "ExamDistribution.Env", row.names = T)

#Numerical Testing####
library(nortest)  # ...from the package 'nortest'...
require(nortest)
ad.test(distribution$density)  # Anderson-Darling test = 0.04
cvm.test(distribution$density) #Cramer-von Mises normality test = 0.03
#Not normal distribution

###Latitudes North and South#### 
south <- data.frame(distribution.ordered[1:9,])
north <- data.frame(distribution.ordered[10:25,])
#Correlation South
cor.test(x=as.numeric(south$latitude), y=south$density, 
         method = "spearman" ) 
plot(x=south$latitude, y=south$density, 
         method = "spearman" ) ##rho= 0.5439, p-value= 0.1301
#Correlation North
cor.test(x=north$latitude, y=north$density, 
    method = "spearman" ) # rho=0.9794118, p-value= <2.2e-16
plot(x=north$latitude, y=north$density)

#Plot both
plot(x=distribution.ordered$latitude,y=distribution.ordered$density)

