working_directory <-'/Users/Daisy/Desktop/MultivariateStats/multivar'

#1.Read Data####
env.density <-read.delim(file = paste0(working_directory, '//env.density.txt'), 
                         h=T, sep ="\t", row.names = 1)
env.density$hemisphere<- ifelse(env.density$Latitude < 0,"south","north")
env.density$hemisphere<- as.factor(env.density$hemisphere)
h.env.density<-env.density[,-3]

#Look the data:
plot(h.env.density$density ~ h.env.density$Temperature)  
plot(h.env.density$density ~ h.env.density$hemisphere)
plot(h.env.density$density ~ h.env.density$Macrophyte)

##2.Tree Model####
library(tree)
tree_model <- tree::tree(density ~ hemisphere+Temperature+Macrophyte,
                         data = h.env.density )
quartz()
plot(tree_model)
text(tree_model)
dev.off()
#3.Ancova model####
model.1 <- lm(formula = density ~ Temperature*Macrophyte*hemisphere, 
            data = h.env.density)
summary.aov(model.1)
summary.lm(model.1)

model.2 <- update(model.1,~.-Temperature:Macrophyte:hemisphere)
car::vif(mod= model.2)
summary.aov(model.2)
summary.lm(model.2)

model.3 <- update(model.2,~.-Temperature:Macrophyte)
car::vif(mod= model.3)
summary.aov(model.3)
summary.lm(model.3)


model.final <- lm(density~Temperature+Macrophyte+hemisphere+
                    Temperature*hemisphere+ Macrophyte*hemisphere, 
                  data = h.env.density)
car::vif(mod= model.final)
summary.aov(model.final)
summary.lm(model.final)

#4. Model checking ####
#Check if the final model is not significantly worse than the maximal model:
AIC(model.1, model.final)     
anova(model.1, model.final)   # not significantly different 

# Check for multicollinearity in this model using VIF.
car::vif(mod = model.final)
quartz()
par(mfrow = c(2,2))
plot(model.2)
dev.off()
par(mfrow = c(1,1)) 
#Are the residuals normally distributed?
ad.test(resid(model.final)) #0.014 
cvm.test(resid(model.final)) #0.0414 
# NOTNormally distributed 
#5. Log Transform density####
h.env.density$log_density <- log10(h.env.density$density)

#Ancova model
model.A <- lm(formula = log_density ~ Temperature*Macrophyte*hemisphere, 
              data = h.env.density)
summary.aov(model.A)
summary.lm(model.A)

model.B <-  update(model.A,~.-Temperature:Macrophyte)
car::vif(mod= model.B)
summary.aov(model.B)
summary.lm(model.B)

model.C <-  update(model.B,~.-Temperature:Macrophyte:hemisphere)
car::vif(mod= model.C)
summary.aov(model.C)
summary.lm(model.C)

model.new<- lm(log_density ~ Temperature+ Macrophyte+hemisphere+
                 Temperature*hemisphere + Macrophyte*hemisphere,
               data=h.env.density)
summary.aov(model.new)
summary.lm(model.new)

#6. Model checking 2.0####
#Check if the final model is not significantly worse than the maximal model:
AIC(model.A, model.new)     
anova(model.A, model.new)   # not significantly different 

ad.test(resid(model.new)) #0.08 
cvm.test(resid(model.new))#0.1117 #Normally distributed

# Check for multicollinearity in this model using VIF.
car::vif(mod = model.new)
quartz()
par(mfrow = c(2,2))
plot(model.C)
dev.off()
par(mfrow = c(1,1)) #everything looks okay 

quartz()

#7.Understanding Model####
#No relationship between denisity and macrophyte cover
plot(h.env.density$log_density ~ h.env.density$Macrophyte, 
     cex = 1, pch = 21, bg = (1 + as.numeric(h.env.density$hemisphere)))
dev.off()
#8.Plot log_Denisty v. Temperature####
plot(h.env.density$log_density ~ h.env.density$Temperature, 
     cex = 1, pch = 21, bg = (1 + as.numeric(h.env.density$hemisphere)))
#North = red,
abline(a = c(model.new$coefficients[1],  # intercept...
             model.new$coefficients[2]),col = 'red')  # ... and slope 
abline(a = c(model.new$coefficients[1]+model.new$coefficients[4],  # intercept...
             model.new$coefficients[2]+model.new$coefficients[5]),col = 'green')

#separately for visualization
quartz()
par(mfrow = c(1,2))

plot(log_density ~ Temperature,data=subset(h.env.density,hemisphere=="north"))
abline(a = c(model.new$coefficients[1],  # intercept...
             model.new$coefficients[2]),col = 'red')

plot(log_density ~ Temperature,data=subset(h.env.density,hemisphere=="south"))
abline(a = c(model.new$coefficients[1]+model.new$coefficients[4],  # intercept...
             model.new$coefficients[2]+model.new$coefficients[5]),col = 'green')


###There is a stronger relationship (more significant) between 
    #sea urchin density and temperature in the northern hemisphere 
    #than in the south 
##Neither hemisphere has a relationship between 
  #sea urchin density and macrophyte cover

##SEE summary.lm(model.new) Line 86



load("Exam C Answer")

