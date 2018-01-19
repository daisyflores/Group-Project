working_directory <-'/Users/Daisy/Desktop/MultivariateStats/multivar'

#Read Data####
env <-read.delim(file = paste0(working_directory, '//data2_env.txt'), 
                 h=T, sep ="\t", row.names = 1)
species <-(read.delim(file = paste0(working_directory, '//data2_species.txt'), 
                      h=T, sep ="\t"))

#Calculate abundance/transect ####
require(reshape) # (load necessary package)
species.wide<-cast(species,siteID+transect~species,sum)
species.melted <- melt(species.wide)
row.names(species.melted) <- 1:length(species.melted[,1])
names(species.melted) <- c(names(species.melted)[1:2],'density',
                           names(species.melted)[4])
species.melted <- subset(species.melted,select=c(1:3))
species.density<- aggregate(density~ siteID,
                            FUN=sum,data=species.melted)

##Prepare to combine into one table ####
latitude<- aggregate(Latitude~ siteID,
                            FUN=mean,env)


all.equal(as.character(latitude$siteID),
          as.character(species.density$siteID))

#Combine Data
distribution <- data.frame(species.density,latitude$Latitude)
distribution.ordered <- distribution[order(distribution$latitude.Latitude),]
row.names(distribution.ordered) <- 1:length(distribution.ordered[,1])

write.table(distribution, file = "Distribution.Env", row.names = T)

#Numerical Testing####
library(nortest)  # ...from the package 'nortest'...
require(nortest)
ad.test(distribution$density)  # Anderson-Darling test and 
cvm.test(distribution$density)

###Latitudes North and South#### 
south <- data.frame(distribution.ordered[1:9,])
north <- data.frame(distribution.ordered[10:25,])
#Correlation South
cor.test(x=as.numeric(south$latitude.Latitude), y=south$density, 
         method = "spearman" ) 
plot(x=south$latitude.Latitude, y=south$density, 
         method = "spearman" ) ##rho= 0.5439, p-value= 0.1301
#Correlation North
cor.test(x=north$latitude.Latitude, y=north$density, 
    method = "spearman" ) # rho=0.9794118, p-value= <2.2e-16
plot(x=north$latitude.Latitude, y=north$density)

#Plot both
plot(x=distribution.ordered$latitude.Latitude,y=distribution.ordered$density)
