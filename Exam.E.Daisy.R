# e) Are sea urchin functional richness and dispersion related to latitude,
# temperature and/or macrophyte cover, and if so, which of these predictors 
# plays a stronger role driving the spatial variability of these indices?
#  HA: Sea urchin community structure is related to temperature and/or 
# macrophyte cover and certain species are most responsible for these pattern.

rm(list=ls()) 

require(reshape)
#1. Compute FD using Multidim FD####

source("quality_funct_space.R")
source("plot_funct_space.R")
source("multidimFD.R")

# Import traits matrix
working_directory <-'/Users/Daisy/Desktop/MultivariateStats/multivar'
traits<- read.table(paste0(working_directory, '//urchingtraits.txt'), 
                    h=T, row.names = 1)
#Delete Species that is not in our abundance data set.
traits <- traits[-8,] #echinotrix.diademis
# Load the dataset
biomass <-(read.delim(file = paste0(working_directory, '//data2_species.txt'), 
                      h=T, sep ="\t"))

# imported the abundance instead of the bites and transform it into matrix
require(reshape)
biomass.wide <- cast(biomass, siteID~species, sum) 
biomass.wide.df <- data.frame(biomass.wide[,-1], row.names = biomass.wide[,1])

abundance <-as.matrix(biomass.wide.df)
#checking if colnames of the species is matching the other df
sum(row.names(traits) %in% colnames(abundance)) == ncol(abundance)
identical(rownames(traits), colnames(abundance))
# computing all functional spaces based on PCoA (up to 6 axes) BECAUSE WE ONLY HAVE 6 TRAITS SEE nbdim=6
# NOTE: you need to be connected to internet for a correct functioning of quality_funct_space function
qual_funct_space<-quality_funct_space(traits, 
                                      traits_weights=NULL, 
                                      nbdim=6, metric="Euclidean",
                                      dendro=F, plot="quality_funct_space")
qual_funct_space$meanSD


# extracting species coordinates in the 4D space
coord_4D<-qual_funct_space$details_funct_space$mat_coord[,1:4]
# saving them
save(coord_4D, file="coord_4D")

# computing functional diversity indices 
FD_baskets<-multidimFD(coord_4D,abundance, check_species_pool=F, verb=TRUE)

# printing results = rounded FD indices values
IndicesAll<-round(FD_baskets,3)
Indices<-as.data.frame(IndicesAll[,c(1,19,21:22)])
Indices #Functional Richness = Column 2 and Functional Dispersion = Column4

# exporting to excel
setwd("/Users/Daisy/Desktop/MultivariateStats/multivar") 
write.table(Indices, file="ObservedIndicesMultidimFD.txt", sep="/t") 



# load the "environmental data":
env.data <- read.table("data2_env.txt", h = T, sep = "\t")

# means of environmental data (temp & macro) by siteide and latitude so we end up with 25 rows
env.data <- aggregate(env.data[, 3:5], list(env.data$siteID), mean)

# rename the first column from Group.1 to site
colnames(env.data)[1] <- "site"

#are the sites in the same order?
Indices <- cbind(site = rownames(Indices), Indices)
identical(env.data$site,Indices$site)

dat <- data.frame(cbind(env.data,Indices[,3:5]))
allbiodiversity<-dat[-1]

# Start of univariate models with Functional Dispersion ####
# (1) Data exploration
# 1.1. Checking if you have outliers in your response variable (Diversity)
par(mfrow = c(1, 2))
boxplot(allbiodiversity$FDis, 
        main = "TaxDiversity")
dotchart(allbiodiversity$FDis, 
         xlab = "Range of data", 
         ylab = "Values")
require(nortest)
# We do not have outliers, these would show as sites 
# far away from others

# 1.2. Checking if you have outliers in your explanatory variables
# lets remind ourselves of the names of variables:
names(allbiodiversity)
par(mfrow = c(1, 3), mar = c(4, 3, 3, 2))
dotchart(allbiodiversity$Latitude, main = "LAT")
dotchart(allbiodiversity$Temperature, main = "TEMP")
dotchart(allbiodiversity$Macrophyte, main = "MACRO")
# temp 1 outlier?

# 1.3. Check for collinearity
panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, use="complete.obs"))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
  
  test <- cor.test(x,y)
  Signif <- symnum(test$p.value, corr=FALSE, na=FALSE, 
                   cutpoints = c(0, 0.05, 0.1, 1), 
                   symbols = c("*", ".", " "))
  
  text(0.5, 0.5, txt, cex = cex * r)
  text (.8, .8, Signif, cex=cex, col=2)
}

# plot every pair of variables, asking for the plot to include correlation coefficients
pairs(allbiodiversity, lower.panel = panel.smooth, upper.panel = panel.cor)

# 1.4. Check for Zero inflation of the response variable (too many zeroes)
sum(allbiodiversity$FDis == 0)
sum(allbiodiversity$FDis == 0) / nrow(allbiodiversity)
plot(table(allbiodiversity$FDis)) #More useful for count data
# not a problem here. not many values are 0 (in the Y axis)

# 2. Fitting the model (including all possible factors and interactions)
# in this case no interactions are possible:
# only consider including interactions when there are sound ecological 
# hypotheses behind them, and when you have enough replication within the
# levels of the interaction. 
modeldiv<-lm(FDis~Temperature+Macrophyte, data=allbiodiversity)
summary(modeldiv) # Latitude effect on edge of significance

# Remove insiginificant term: Macrophyte
modeldiv2 <- update(object = modeldiv,  ~ . - Macrophyte)
summary(modeldiv2) 

# Check residuals:
par(mfrow = c(2,2))
plot(modeldiv2)
# Are the residuals normally distributed?

# The model should show homogeneity of variance (heteroscedasticity):
# (i.e. the first plot residuals vs. fitted should show no pattern
# all dots spread across all biplot)
par(mfrow = c(2, 2))
plot(modeldiv)
resid.vals.mod <- resid(modeldiv)
fitted.vals.mod <- fitted(modeldiv)
plot(fitted.vals.mod,resid.vals.mod) # not sure what goes in here, not clear 
abline(h = 0, v = 0, lty = 2)         # from sonias script just said x = F2, y = E2
# all dots are spread across plot, is homogneous

#WHAT DO WE CONCLUDE????
# Richness ####
# 1.1. Checking if you have outliers in your response variable
par(mfrow = c(1, 2))
boxplot(allbiodiversity$FRic, 
        main = "Functional Richness")
dotchart(allbiodiversity$FRic, 
         xlab = "Range of data", 
         ylab = "Values")
# No Outliers


# 1.4. Check for Zero inflation of the response variable (too many zeroes)
sum(allbiodiversity$FRic == 0)
sum(allbiodiversity$FRic == 0) / nrow(allbiodiversity)
plot(table(allbiodiversity$FRic))
# no 0 inflation

# 2. Fitting the model (including all possible factors and interactions)
modeldiv<-lm(FRic~Temperature+Macrophyte, data=allbiodiversity)
summary(modeldiv) 
# all are insignificant

# as expected the model has a very bad Rsquared, and is not significant, 
# but before we trust it lets run diagnostic plots:

par(mfrow = c(2, 2))
plot(modeldiv)
resid.vals.mod <- resid(modeldiv)
fitted.vals.mod <- fitted(modeldiv)
plot(x = fitted.vals.mod, y = resid.vals.mod)
abline(h = 0, v = 0, lty = 2)
plot(modeldiv)
# all dots are spread across plot, is homogneous

# functional richness and dispersion are colinear 
#They are NOT related to any of the environmental variables 
