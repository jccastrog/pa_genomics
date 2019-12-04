#!/usr/bin/env Rscript
#plotTraits.R

args = commandArgs(TRUE)
effects <- read.table(args[1], h = T, stringsAsFactors = F)
effects <- read.table("PA_trait_4_effects.tsv", h = T, stringsAsFactors = F)

#Manhattan plot
plot(-log(effects$mean.p.value, base = 10), col = 'grey55', xlab='', ylab = expression('-log'[10]*'(P)'))
#Volcano plot
plot(effects$mean.effects,-log(effects$mean.p.value,base = 10),pch=19,col='grey55',xlab='Effect',ylab=expression('-log'[10]*'(P)'),main='Biofilm')


for (i in 3:16){
  #file <- paste("PA_trait_", i, "_effects.tsv", sep = "")
  #effects <- read.table(file, h = T, stringsAsFactors = F)
  #Volcano plot
  #plot(effects$mean.effects,-log(effects$mean.p.value,base = 10),pch=19,col='grey55',xlab='Effect',ylab=expression('-log'[10]*'(P)'),main='Biofilm')
  file <- paste("PA_trait_", i, "_mse.tsv", sep = "")
  mse <- read.table(file, h = T, stringsAsFactors = F)
  #Boxplot MSE
  boxplot(mse)
}

y <- phenotype.data$biofilm.live
X <- clean.data
saveAt <- tempdir()
df0 <- 5
S0 <- NULL
weights <- NULL

ETA <- list(list(X = X, model = "BayesB"), list(X = X, model = "BL"))
log <- captureOutput({fit_BGLR <- BGLR(y=y,ETA=ETA,nIter=num.iter,burnIn=burn.in,thin=thin,saveAt=saveAt,df0=5,S0=S0,weights=weights,R2=R2)})
rem.str <- paste('rm ',saveAt,'*.dat',sep = '')
system(rem.str)

eff <- fitBGLR$ETA[[1]]$b
pVal <- pt(-abs(eff/fitBGLR$ETA[[1]]$SD.b), df = 5)
plot(eff,-log(pVal,base = 10),pch=19,col='grey55',xlab='Effect',ylab=expression('-log'[10]*'(P)'),main=alg)
