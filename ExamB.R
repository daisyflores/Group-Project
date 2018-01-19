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
species.density<- aggregate(abundance~ siteID+transect,
                            FUN=sum,data=species.melted)
species.density$density <- species.density$abundance/50


##3.Prepare to combine into one table ####
density.ordered <- species.density[
  order(species.density$siteID, species.density$transect),
  ]

env.ordered <- env[order(env$siteID, env$transect),
  ]

#Determine whether they are ordered the same before combining data 
all.equal(as.character(density.ordered$siteID),
          as.character(env.ordered$siteID)) #True
#Combine Data
env.density <- data.frame(env.ordered,density.ordered$density)
names(env.density)[6] <- c("density")
env.density <- env.density[order(env.density$Latitude),]
row.names(env.density) <- 1:length(env.density[,1])
write.table(env.density, 
            "/Users/Daisy/Desktop/MultivariateStats/multivar/env.density.txt", 
            sep="\t")
#3b. Absolute Value?####
env.density$Latitude <- as.numeric(lapply(env.density$Latitude, abs))

#4.Look at all correlations####
quartz()
#pairs(env.density,panel=panel.smooth)
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
pairs(env.density[ ,3:6], lower.panel = panel.smooth, upper.panel = panel.cor)
dev.off()
#HELP GAM MODEL ####

gam_model <- gam(formula = density ~ s(Temperature)+s(Macrophyte),
                 data = env.density,
                 family = 'gaussian',
                 link = 'identity')


#  Now plot the GAM model:
quartz()
par(mfrow = c(2,2))
plot(gam_model)
dev.off()
par(mfrow=c(1,1))

#Confidence intervals are sufficiently narrow= curvature of density response


#5.Fit a tree-model#####
#to search for complex interactions between expl vars
library(tree)
tree_model <- tree::tree(density ~ Temperature+Macrophyte,
                         data = env.density)
quartz()
plot(tree_model)
text(tree_model)
dev.off()
#Complicated relationship between Temperature and Macrophyte: gonna be hard modellll

    


#7. Construct model, with all factors, interactions and covariates ####
#       Interactions between all three explanatory variables,
#       plus quadratic terms to test for curvature in response to each of the three
#       explanatory variables.
model.1<- lm(density~Temperature*Macrophyte + 
              I(Temperature^2)+ I(Macrophyte^2), 
            data = env.density)
require(car)
car::vif(mod=model.1)

summary(model.1)

#Simplify model####
      #Step term####
model.2a<-step(model.1)
summary(model.2a)

#Manually simplified#
model.2<- update(model.1,~.-Temperature:Macrophyte)
car::vif(mod= model.2)
summary(model.2)
model.3<- update(model.2,~.-I(Macrophyte^2))
car::vif(mod= model.3)
summary(model.3)


#8. Model Checking####
#Check if the final model is not significantly worse than the maximal model:
AIC(model.1, model.2a)     
anova(model.1, model.2a)   # not significantly different 

# Check for multicollinearity in this model using VIF.
car::vif(mod = model.2a)
quartz()
par(mfrow = c(2,2))
plot(model.2a)
dev.off()
par(mfrow = c(1,1)) 
#Are the residuals normally distributed now?
ad.test(resid(model.2a)) #0.0667
cvm.test(resid(model.2a)) #0.2287
# Normally distributed 



##For Report:
      #ANOVA tables are often published containing a
      #mixture of significant and non-significant effects.
      #say whether your data are orthogonal or not
      #present a minimal adequate model
      #give a list of the non-significant terms that were omitted, 
      #and the deviance changes that #resulted from their deletion

