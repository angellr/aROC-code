
################################
##
## Written by Bob Angell, University of Utah School of Medicine, Biomedical Informatics, 2018
## This code may be freely used if proper attribution is used. 
## 
## Contact Information: bob.angell@utah.edu; angellr@mac.com; mobile: 801-706-2520
##
################################

################################
##
## This code will create a single graph plotting three ROC (AUC) curves. Below, one needs to read
## in the three datasets for this analysis. I have provided (3) anonymized datasets as examples.
##
################################

library(pROC)
library(Hmisc)
library(data.table)
library(ROCR)
library(grDevices)

ResultsOfRun1 <-  Neuropathy
setnames(ResultsOfRun1, "Neuropathy", "ReferenceState")
setnames(ResultsOfRun1, "ProbNeuropathy", "Probability")
describe(ResultsOfRun1)

ResultsOfRun2 <-  Retinopathy
setnames(ResultsOfRun2, "Retinopathy", "ReferenceState")
setnames(ResultsOfRun2, "ProbRetinopathy", "Probability")
describe(ResultsOfRun2)

ResultsOfRun3 <-  Nephropathy
setnames(ResultsOfRun3, "Nephropathy", "ReferenceState")
setnames(ResultsOfRun3, "ProbNephropathy", "Probability")
describe(ResultsOfRun3)

###############################
##
## Now, lets create the 3 ROC models
##
###############################

roc.model1 = roc(ResultsOfRun1$ReferenceState, ResultsOfRun1$Probability, levels=c("Absent", "Present"))
roc.model2 = roc(ResultsOfRun2$ReferenceState, ResultsOfRun2$Probability, levels=c("Absent", "Present"))
roc.model3 = roc(ResultsOfRun3$ReferenceState, ResultsOfRun3$Probability, levels=c("Absent", "Present"))

###############################
##
## Calculate the area under the receiver operating curve for the 3 models.
##
###############################

auc(roc.model1)
auc(roc.model2)
auc(roc.model3)

###############################
##
## Calculate the 95% confidence intervals for the 3 models.
##
###############################

ci.auc(roc.model1)
ci.auc(roc.model2)
ci.auc(roc.model3)

#######################################################################################
## Various Plotting Options
## Handy if test set has small numbers.
#######################################################################################

plot.roc(roc.model1, col = "black") ## no smoothing
lines(smooth(roc.model1, method="density"), col = "purple") ## Density
lines(smooth(roc.model1), method="binormal", col = "red") ## Binormal smoothing

plot.roc(roc.model2, col = "black") ## no smoothing
lines(smooth(roc.model2, method="density"), col = "purple") ## Density
lines(smooth(roc.model2), method="binormal", col = "red") ## Binormal smoothing

plot.roc(roc.model3, col = "black", add=TRUE) ## no smoothing
lines(smooth(roc.model3, method="density"), col = "purple") ## Density
lines(smooth(roc.model3), method="binormal", col = "red") ## Binormal smoothing


##############################################
##
## Regular Plotting Options
## ROC and 95% Confidence Levels.
##
##############################################

##Plot simple ROC model to get a look
plot.roc(roc.model1, col = "black", print.auc=TRUE, grid=c(0.1, 0.1), grid.col=c("green", "red"))

plot.roc(roc.model2, col = "black", print.auc=TRUE, grid=c(0.1, 0.1), grid.col=c("green", "red"))

plot.roc(roc.model3, col = "black", print.auc=TRUE, grid=c(0.1, 0.1), grid.col=c("green", "red"))


##find AUC for model and find confidence interval for AUC
auc(roc.model1)
ci.auc(roc.model1)

auc(roc.model2)
ci.auc(roc.model2)

auc(roc.model3)
ci.auc(roc.model3)

###############################
##
## This will graph each model and draw 95% CI around the curve.
##
###############################

##find and plot confidence intervals
shape.cis1 <- ci.sp(roc.model1, sensitivities=seq(0, 1, .01), boot.n=500)
print(shape.cis1)
plot(shape.cis1, type="shape", col="yellow")

shape.cis2 <- ci.sp(roc.model2, sensitivities=seq(0, 1, .01), boot.n=500)
print(shape.cis2)
plot(shape.cis2, type="shape", col="pink")

shape.cis3 <- ci.sp(roc.model3, sensitivities=seq(0, 1, .01), boot.n=500)
print(shape.cis3)
plot(shape.cis3, type="shape", col="salmon")

###############################
##
## Table of Confidence Intervals
## ROC and 95% Confidence Levels.
##
###############################

ci(roc.model1)
auc(roc.model1)
thresholds <- ci.thresholds(roc.model1, boot.stratified=TRUE, conf.level=0.95, boot.n=2000, thresholds=seq(0, 1, 0.01))
print(thresholds)

ci(roc.model2)
auc(roc.model2)
thresholds <- ci.thresholds(roc.model2, boot.stratified=TRUE, conf.level=0.95, boot.n=2000, thresholds=seq(0, 1, 0.01))
print(thresholds)

ci(roc.model3)
auc(roc.model3)
thresholds <- ci.thresholds(roc.model3, boot.stratified=TRUE, conf.level=0.95, boot.n=2000, thresholds=seq(0, 1, 0.01))
print(thresholds)

###############################

###############################
##
## Set the colors we want for the 3 models.
##
###############################

color1 <- adjustcolor("black", alpha.f = 1.0)
color2 <- adjustcolor("red", alpha.f = 1.0)
color3 <- adjustcolor("blue", alpha.f = 1.0)

color4 <- adjustcolor("lightsalmon", alpha.f = 0.8)
color5 <- adjustcolor("pink", alpha.f = 0.8)
color6 <- adjustcolor("yellow", alpha.f = 0.8)

color7 <- adjustcolor("red", alpha.f = 0.6)
color8 <- adjustcolor("green", alpha.f = 0.6)

###############################
##
## Here is where the magic happens with 3-in-1 ... this is a nice graph.
##
###############################

plot(smooth(roc.model1), col = color1, print.auc=TRUE, print.auc.y = .6, lwd=3, main = "aROC Microvasculature Model Results",grid=c(0.1, 0.1), grid.col=c(color7, color8))
plot(smooth(roc.model2), col = color2, print.auc=TRUE, print.auc.y = .5, lwd=3, add = TRUE)
plot(smooth(roc.model3), col = color3, print.auc=TRUE, print.auc.y = .4, lwd=3, add = TRUE)
plot(legend("bottomright", legend=c("Neuropathy", "Retinopathy", "Nephropathy"),col=c(par("fg"), color2, color3, color1), lwd=3), add=TRUE)


