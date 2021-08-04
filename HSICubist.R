#install.packages("caret")
#install.packages("Cubist")
#install.packages("chillR")
#install.packages("openxlsx")
#install.packages("readxl")
#install.packages("lattice")
#install.packages("ithir")

library(prospectr)
library(pls)
library(chillR)
library(caret)
library(Cubist)
library(ggplot2)
library(readxl)
library(openxlsx)
library(ithir)
library (kernlab)

cat("\f") 
rm(list=ls())

############## setting up the data set ############## 
setwd('C:\\Users\\pascalpc.com\\OneDrive\\Czech\\Asa\\VegPTE\\Data\\Codes') 
dat<- read.csv(file='AllData.csv', sep=",",header=T, row.names=1)
str(dat)

dat = dat[complete.cases(dat),]   ##removes the rows with NA values in any of the columns
spectra1<- dat[, c(7:ncol(dat))]  ##subset the data set which contains only the spectra (raw spectra)
spectra<- dat[, c(16:763)]        ##subset the data set which contains only the selected region (from 390 to 1090)
names(spectra) <- gsub("X", "", names(spectra)) ## Across all columns, replace all instances of "X" with ""
str(spectra)  ##prints the structure of the subset data set

HSI.spectra<- as.matrix(spectra)
dat.Lab <- dat[,c(1:6)]  ##Reference data set (lab measurements)
dat.all<-data.frame(dat.Lab, HSI.spectra)
wav <- as.numeric(colnames(spectra))
matplot(wav, t(HSI.spectra), type="l", lty=1, xlab="Wavelength (nm)", ylab= "Absorbance (-)", col="blue")



############## PCA ############## 
HSI.spectra.pca <- prcomp(HSI.spectra,scale =FALSE)
str(HSI.spectra.pca)
plot(HSI.spectra.pca)

eig <- (HSI.spectra.pca$sdev)^2  ##Eigenvalues
variance <- eig*100/sum(eig) ##Variances in percentage
cumvar <- cumsum(variance) ##Cumulative variances
dat.var <- data.frame(eig = eig, variance = variance,
                      cumvariance = cumvar)
head(dat.var)

### score plot
plot(c(1:10), dat.var[1:10, 2],
     main = "Variances",
     xlab = "Principal Components (PCs)",
     ylab = "Percentage of variances",
     xaxs = 'i',
     yaxs = 'i',
     xlim = c(0,10),
     ylim = c(0,100))
### Add connected line segments to the plot
lines(x = 1:10, dat.var[1:10, 2], col = "red")

### PCA score plot 
plot(HSI.spectra.pca$x[,1:2], pch=19, xlab = paste('PC1 (', round(dat.var$variance[1],2),'%)', sep=''), 
     ylab = paste('PC1 (', round(dat.var$variance[2],2),'%)', sep=''))

abline(h=0,v=0, lty=2, col="gray")

### PCA loading plot
matplot(wav, -HSI.spectra.pca$rotation[,1:2], type ="l",lty= 1:2, col = 1:2,xlab = "Wavelengths (nm)", ylab = "Loadings (-)")
abline(h=0, col="gray", lty=2)
legend("topright", legend = paste(c("PC1", "PC2"), c("(","("),round(variance[1:2],2),c("%)", "%)")) , col=1:2, lty=1:2)



############## Cook Distance ############## 
### Finding outliers
mod <- lm(dat.Lab, data=dat.all)
cooksd <- cooks.distance(mod)

### Plot the Cook's Distance using the traditional 4/n criterion
sample_size <- nrow(HSI.spectra)
plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance")  ##plot cook's distance
abline(h = 4/sample_size, col="red")  ##add cutoff line
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4/sample_size, names(cooksd),""), col="red")  ##add labels

### Removing outliers ########
### influential row numbers
influential <- as.numeric(names(cooksd)[(cooksd > (4/sample_size))])
Outliers <- HSI.spectra[-influential, ]

plot3 <- ggplot(HSI.spectra, aes(x = speed, y = dist)) +
  geom_point() + 
  geom_smooth(method = lm) +
  xlim(0, 20) + ylim(0, 220) + 
  ggtitle("Before")
plot4 <- ggplot(Outliers, aes(x = speed, y = dist)) +
  geom_point() + 
  geom_smooth(method = lm) +
  xlim(0, 20) + ylim(0, 220) + 
  ggtitle("After")

gridExtra::grid.arrange(plot3, plot4, ncol=2)


############## Pre-processing ############## 
sg <- savitzkyGolay(HSI.spectra,p=2,w=11,m=2) ##savitzky-golay:Run this regardless of which other preprocessing technique are applied 
#der <- t(diff(t(HSI.spectra), differences = 1, lag = 11)) ##derivative
#gder<- gapDer(HSI.spectra, m = 2, w = 11, s = 2)##gap derivative
#cr <- continuumRemoval(HSI.spectra,wav,type='A') ##continuum-removal
#snv <- standardNormalVariate(HSI.spectra)   ##Standard Normal Variate 
#msc <- msc(HSI.spectra)   ##Multiplicative Scatter Correction
msc <- msc(sg)   ##based on sg

dat.n<-data.frame(dat.Lab, msc=I(msc))
matplot(wav, t(dat.n$msc), type="l", lty=1, xlab="Wavelength (nm)", ylab= "Absorbance (-)", col="red")


############ Modeling ############
####Split data (Kennard Stone)
#ken.sto(inp, per = "T", per.n = 0.3, num, va = "F", sav = "T", path = "", out = "Sel")
cDat<- kenStone(HSI.spectra, per = "T", per.n = 0.7, num, va = "F", sav = "T", path = "C:\\Users\\pascalpc.com\\OneDrive\\Czech\\Asa\\VegPTE\\Data\\Codes", out = "KS_splitSamples")


### Fit the model
### Cubist #change the input based on the KS output
nc<-ncol(cDat) 
ctrl<- trainControl(method = "cv", number= 5, #K-fold cross-validation
                    summaryFunction = defaultSummary, selectionFunction = "best")

fit_cubist<-train(x = cDat[, 7:nc],
                  y = cDat$Fe,method="cubist",trControl=ctrl)
summary(fit_cubist)
print(fit_cubist)
### Cubist validation (R2, RMSE, MSEM bias)
Cubist.pred.V<-predict(fit_cubist,newdata=vDat)
goof(observed = vDat$Fe,predicted = Cubist.pred.V)

### SVM
set.seed(1)
fit_SVM <- train(x = cDat[, 7:nc], y = cDat$Fe, method = "svmRadial", trContro = ctrl)
summary(fit_SVM)
print(fit_SVM)
### SVM validation (R2, RMSE, MSEM bias)
SVM.pred.V<-predict(fit_cubist,newdata=vDat)
goof(observed = vDat$Fe,predicted = SVM.pred.V)               
