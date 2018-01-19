rm(list=ls()) 
require(reshape) 
require(ade4)
require(vegan)
require(packfor)
require(MASS)
require(ellipse)
require(FactoMineR)

# Set the working directory:
setwd("/Users/Daisy/Desktop/MultivariateStats/multivar")

biomass <- read.table("data2_species.txt", h= T, sep = "\t")
biomass.wide<-cast(biomass,siteID~species,sum)
biomass.wide.df <- data.frame(biomass.wide[,-1], row.names=biomass.wide[,1])

env <- read.table("data2_env.txt", h= T,sep = "\t")
env.mean<- aggregate(list(Latitude=env$Latitude, Temperature=env$Temperature,
                          Macrophyte=env$Macrophyte),
                     by = list(siteID=env$siteID),mean)
env.mean <- data.frame(env.mean[,-1], row.names=env.mean[,1])
subset.env.mean<-env.mean[c(F, T, T)] 

# Load some graphic functions designed by Brocard et al (2011) we will need:
source("evplot.R")
source("cleanplot.pca.R")

# Hellinger-transform the species dataset
spe.hel<-decostand(biomass.wide.df, "hellinger")

# standardise the reduced environmental dataset:
stand.env.mean<-scale(subset.env.mean)
stand.env.mean<-as.data.frame(stand.env.mean)

# Run the simple RDA
formulaRDA<-rda(spe.hel ~. , data = stand.env.mean)
summary(formulaRDA)

#Significance Test
anova(formulaRDA)
          #Df Variance      F     Pr(>F)    
# Model     2  0.10070    19.756  0.001 ***
# Residual 22  0.05607

#It is also possible to analyse terms separately:
anova(formulaRDA, by="term", permu=600) 
#Or analyse significances of marginal effects (“Type III effects”)
anova(formulaRDA, by="mar")
#Analyse significance of each axis
anova(formulaRDA, by="axis", perm=500)

# extract the canonical coefficient from RDA object:
coef(formulaRDA)


# Check whether the whole model is statistically significant:
set.seed(111)
anova.cca(formulaRDA, step=1000)
# is significant p=0.001

# compute the R2 and adjusted R2
(R2<-RsquareAdj(formulaRDA)$r.squared) # 0.6423466
(R2adj <-RsquareAdj(formulaRDA)$adj.r.squared) #0.6098326
# The second value is the unbiased amount of variation explained 
# by the model

# Test the significance of all the canonical axes
set.seed(111)
anova.cca(formulaRDA, by = "axis", step=1000)


quartz()
#from sonia: Plot1 (this scaling should be used to focus mainly the angles of environmental vectors)
plot(formulaRDA, scaling=1, main="RDA analysis of ...")
spe.sc<-scores(formulaRDA, choices=1:2, scaling=1, display="sp")
arrows(0, 0, spe.sc[,1], spe.sc[,2], length = 0, lty=1, col="red")

dev.off()
quartz()
plot(formulaRDA, scaling=1, display=c("sp", "lc", "cn"),
     main="triplot with LC scores")
arrows(0, 0, spe.sc[, 1], spe.sc[ ,2], length=0, lty=1, col="red")
