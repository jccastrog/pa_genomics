suppressPackageStartupMessages(library(BGLR))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(lars))

rm(list = ls())
#load('F:/JuanCamilo/GaTech/BrownLab/projects/pa_genomics/genotype/myEnvironment.RData')
load('/media/jccastrog/Data/JuanCamilo/GaTech/BrownLab/projects/pa_genomics/genotype/myEnvironment.RData')

i <- 1
trait <- 1
mse.table <- data.frame(fold = c(), test = c(), all = c())
params.effects <- data.frame(row.names = og.names)
params.sd <- data.frame(row.names = og.names)
pheno.fold <- phenotype.data
pheno.fold[folds.matrix[i,],] <- NA
y.fold <- as.numeric(pheno.fold[[trait]])
y.test <- as.numeric(phenotype.data[[trait]])

mse.vec <- c()
mse.vec.test <- c()
mse.min <- 400
mse.max <- 0
mse.cv <- 2000
for (q in 1:100){
  fit.BGLR <- runBGLR(X.fold, y.fold, alg = "BL")
  mse.fold <- mean((fit.BGLR$y - fit.BGLR$yHat)^2, na.rm = TRUE)
  mse.test <- mean((y.test[folds.matrix[i,]] - predict(fit.BGLR)[folds.matrix[i,]])^2)
  mse.vec <- c(mse.vec, mse.fold)
  mse.vec.test <- c(mse.vec.test, mse.test)
  if (min(mse.vec) < mse.min){
    mse.min <- min(mse.vec)
    min.fit <- fit.BGLR
  }
  if (max(mse.vec) > mse.max){
    mse.max <- max(mse.vec)
    max.fit <- fit.BGLR
  }
  if (min(mse.vec.test) < mse.cv){
    mse.cv <- min(mse.vec.test)
    min.cv <- fit.BGLR
  }
}
png('/media/jccastrog/Data/JuanCamilo/GaTech/BrownLab/projects/pa_genomics/plots/histogram_MSEbiofilm.pdf')
hist(mse.vec, xlab = 'MSE', col = 'grey55', main ='Fold MSE', breaks = 50)
hist(mse.vec.test, xlab = 'MSE', col = 'grey55', main ='Test MSE', breaks = 50)
graphics.off()

pdf('/media/jccastrog/Data/JuanCamilo/GaTech/BrownLab/projects/pa_genomics/plots/plot_PredictBiofilm.pdf')
plot(min.fit$y, min.fit$yHat, pch = 19, col = 'green', xlab = 'y', ylab = expression(hat(y)))
points(y.test[folds.matrix[i,]] , predict(min.fit)[folds.matrix[i,]], pch = 19, col = 'red')

plot(min.cv$y, min.cv$yHat, pch = 19, col = 'green', xZzlab = 'y', ylab = expression(hat(y)))
points(y.test[folds.matrix[i,]] , predict(min.cv)[folds.matrix[i,]], pch = 19, col = 'red')


plot(max.fit$y, max.fit$yHat, pch = 19, col = 'green', xlab = 'y', ylab = expression(hat(y)))
points(y.test[folds.matrix[i,]] , predict(max.fit)[folds.matrix[i,]], pch = 19, col = 'red')
graphics.off()