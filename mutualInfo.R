#!/usr/bin/env Rscript

################################################################################
# Name:	    mutualInfo.R
# Author:   Juan C. Castro [jcastro37@gatech.edu]
#           Georgia Institute of Technology
#           
# Version:  1.0
# Date:     04-Jan-2018
# License:  GNU General Public License v3.0.
# ==============================================================================
# 
################################################################################

#======================= 0.0 Install required packages ========================#
personal.lib.path = Sys.getenv("R_LIBS_USER")
if(!file.exists(personal.lib.path))
  dir.create(personal.lib.path)

packages <- c("optparse", "entropy", "igraph")
if(any(!(packages %in% installed.packages()))){
  cat("Please wait a moment! Installing required packages ...\n")
  install.packages(packages[!(packages %in% installed.packages())],
                   quiet = T, repos="http://cran.rstudio.com/",
                   lib = personal.lib.path)
  cat("Required packages installed!\n")
}
#======== 1.0 Load packages, define functions, and initialize variable========#
# 1.1 Load packages ==========================================================#
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(entropy))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(lars))
# 1.2 Define functions =======================================================#
#' Estimate mutual information (MI) for all pairs of variables in an expression
#' matrix
#'
#' @param matrixData A matrix with expression data where columns are genes and
#'        rows are time points
#' @param numGenes The number of genes in the matrix
#' @return mutualMat A matrix with values of mutual information for all pair of
#'        genes 
#' @author Juan C. Castro \email{jcastro37@gatech.edu}
mutualInfoEst <- function(matrixData,numGenes) {
  numBins <- ceiling(log2(nrow(matrixData)+1))
  mutualMat <- matrix(ncol=numGenes,nrow=numGenes)
  for (i in 1:numGenes) {
    for (j in 1:numGenes) {
      numVarsI <- length(table(matrixData[,i]))
      numVarsJ <- length(table(matrixData[,j]))
      if (numVarsI >1 & numVarsJ >1){
        discretVec <- discretize2d(matrixData[,i],matrixData[,j],numBins,numBins)
        mutualMat[i,j] <- suppressWarnings(mi.empirical(discretVec,unit=c("log2")))
      } else {
        mutualMat[i,j] <- 0
      }
    }
  }
  return(mutualMat)
}
#' Estimate a null distribution of MI for variables in an expression matrix
#'
#' @param matrixData A matrix with expression data where columns are genes and
#'        rows are time points
#' @param numGenes The number of genes in the matrix
#' @param randomizations Number of times to randomize expression values
#' @return distMat A matrix with values of mutual information for a randomized
#'        expression matrix where each column is a linearized version of a random
#'        expression matrix
#' @author Juan C. Castro \email{jcastro37@gatech.edu}
nullInfoDist <- function(matrixData,numGenes,randomizations){
  distMat <- matrix(ncol=(numGenes^2),nrow=randomizations)
  for (counts in 1:randomizations) {
    randomMatrixData <- matrix(ncol=ncol(matrixData),nrow=(nrow(matrixData)))
    for (i in 1:ncol(matrixData)){
      randomMatrixData[,i] <- sample(matrixData[,i])
    }
    cat(paste("Randomization ",counts,"\n", sep =""))
    randomMIMat <- mutualInfoEst(randomMatrixData,numGenes)
    randomMIVec <- matrix(randomMIMat,ncol=(numGenes^2))
    distMat[counts,] <- randomMIVec
  }
  return(distMat)
}
#' Calculate a score for each value of mutual information in an MI matrix
#'
#' @param initialMatrix A matrix with MI values calculated form an expression
#'        matrix
#' @param nullDistMaxtrix A matrix with null values of MI obatined with nullInfoDist
#' @param numGenes The number of genes in the matrix
#' @param randomizations Number of times to randomize expression values
#' @return pMat A matrix with pValues of for the MI values in initialMatrix
#' @author Juan C. Castro \email{jcastro37@gatech.edu}
infoScore <- function(initialMatrix,nullDistMatrix,numGenes,randomizations){
  matSize <- numGenes^2
  sumDist <- randomizations
  linealIniMat <- matrix(initialMatrix,ncol=matSize)
  pVec <- c()
  for (i in 1:matSize){
    MICounts <- sum(nullDistMatrix[,i] >= linealIniMat[i])
    pScore <- MICounts/sumDist
    pVec[i] <- pScore
  }
  pMat <- matrix(pVec,ncol=numGenes,nrow=numGenes)
  return(pMat)
}

pangenome.df <- read.table('~/projects/pa_genomics/genotype/pa_genomes.ogs.codes.tsv', 
                           h = T , stringsAsFactors = F, sep = '\t', row.names = 1)
p.vals <- read.table('~/projects/pa_genomics/phylogeny_corrs/pvals.out', h = T, stringsAsFactors = F, sep = '\t', row.names = 1)

# Subset the p values
p.vals.biofilm.live <- subset(p.vals, p.vals$biofilm.live <= 0.001)
p.vals.biofilm.dead <- subset(p.vals, p.vals$biofilm.dead <= 0.001)
p.vals.biofilm.od <- subset(p.vals, p.vals$biofilm.od <= 0.001)
p.vals.biofilm.eps <- subset(p.vals, p.vals$biofilm.eps <= 0.001)
p.vals.fraction.biofilm <- subset(p.vals, p.vals$fraction.biofilm <= 0.001)
p.vals.max.od <- subset(p.vals, p.vals$max.od <= 0.001)
p.vals.max.rate <- subset(p.vals, p.vals$max.rate <= 0.001)
p.vals.elastase <- subset(p.vals, p.vals$elastase <= 0.001)
p.vals.carbenicillin.ratio <- subset(p.vals, p.vals$carbenicillin.ratio <= 0.001)
p.vals.tobramycin.ratio <- subset(p.vals, p.vals$tobramycin.ratio <= 0.001)

# Subset the OGs matrix
biofilm.live.matrix <- subset(pangenome.df, rownames(pangenome.df)%in%row.names(p.vals.biofilm.live))
biofilm.dead.matrix <- subset(pangenome.df, rownames(pangenome.df)%in%row.names(p.vals.biofilm %>% .dead))
biofilm.od.matrix <- subset(pangenome.df, rownames(pangenome.df)%in%row.names(p.vals.biofilm.od))
biofilm.eps.matrix <- subset(pangenome.df, rownames(pangenome.df)%in%row.names(p.vals.biofilm.eps))
fraction.biofilm.matrix <- subset(pangenome.df, rownames(pangenome.df)%in%row.names(p.vals.fraction.biofilm))
max.od.matrix <- subset(pangenome.df, rownames(pangenome.df)%in%row.names(p.vals.max.od))
max.rate.matrix <- subset(pangenome.df, rownames(pangenome.df)%in%row.names(p.vals.max.rate))
elastase.matrix <- subset(pangenome.df, rownames(pangenome.df)%in%row.names(p.vals.elastase))
carbenicillin.matrix <- subset(pangenome.df, rownames(pangenome.df)%in%row.names(p.vals.carbenicillin.ratio))
tobramycin.matrix <- subset(pangenome.df, rownames(pangenome.df)%in%row.names(p.vals.tobramycin.ratio))



# Initial MI matrices
biofilm.live.MI <- mutualInfoEst(matrixData = t(biofilm.live.matrix), numGenes = nrow(biofilm.live.matrix))
biofilm.dead.MI <- mutualInfoEst(matrixData = t(biofilm.dead.matrix), numGenes = nrow(biofilm.dead.matrix))
biofilm.od.MI <- mutualInfoEst(matrixData = t(biofilm.od.matrix), numGenes = nrow(biofilm.od.matrix))
biofilm.eps.MI <- mutualInfoEst(matrixData = t(biofilm.eps.matrix), numGenes = nrow(biofilm.eps.matrix))
fraction.biofilm.MI <- mutualInfoEst(matrixData = t(fraction.biofilm.matrix), numGenes = nrow(fraction.biofilm.matrix))
max.od.MI <- mutualInfoEst(matrixData = t(max.od.matrix), numGenes = nrow(max.od.matrix))
max.rate.MI <- mutualInfoEst(matrixData = t(max.rate.matrix), numGenes = nrow(max.rate.matrix))
elastase.MI <- mutualInfoEst(matrixData = t(elastase.matrix), numGenes = nrow(elastase.matrix))
carbenicillin.MI <- mutualInfoEst(matrixData = t(carbenicillin.matrix), numGenes = nrow(carbenicillin.matrix))
tobramycin.MI <- mutualInfoEst(matrixData = t(tobramycin.matrix), numGenes = nrow(tobramycin.matrix))

# Null MI distributions
biofilm.live.null <- nullInfoDist(matrixData = t(biofilm.live.matrix), numGenes =  nrow(biofilm.live.matrix), randomizations = 1000)
biofilm.dead.null <- nullInfoDist(matrixData = t(biofilm.dead.matrix), numGenes =  nrow(biofilm.dead.matrix), randomizations = 1000)
biofilm.od.null <- nullInfoDist(matrixData = t(biofilm.od.matrix), numGenes =  nrow(biofilm.od.matrix), randomizations = 1000)
biofilm.eps.null <- nullInfoDist(matrixData = t(biofilm.eps.matrix), numGenes =  nrow(biofilm.eps.matrix), randomizations = 1000)
fraction.biofilm.null <- nullInfoDist(matrixData = t(fraction.biofilm.matrix), numGenes =  nrow(fraction.biofilm.matrix), randomizations = 1000)
max.od.null <- nullInfoDist(matrixData = t(max.od.matrix), numGenes =  nrow(max.od.matrix), randomizations = 1000)
max.rate.null <- nullInfoDist(matrixData = t(max.rate.matrix), numGenes =  nrow(max.rate.matrix), randomizations = 1000)
carbenicillin..null <- nullInfoDist(matrixData = t(carbenicillin.matrix), numGenes =  nrow(carbenicillin.matrix), randomizations = 1000)
tobramycin.null<- nullInfoDist(matrixData = t(tobramycin.matrix), numGenes =  nrow(tobramycin.matrix), randomizations = 1000)

# P values for MI
p.values < - infoScore(biofilm.live.MI, biofilm.live.null, nrow(biofilm.live.matrix),1000)
p.values <- infoScore(biofilm.dead.MI, biofilm.dead.null, nrow(biofilm.dead.matrix),1000)
p.values <- infoScore(biofilm.od.MI, biofilm.od.null, nrow(biofilm.od.matrix),1000)
p.values <- infoScore(biofilm.eps.MI, biofilm.eps.null, nrow(biofilm.eps.matrix),1000)
p.values <- infoScore(fraction.biofilm.MI, fraction.biofilm.null, nrow(fraction.biofilm.matrix),1000)
p.values <- infoScore(max.od.MI, max.od.null, nrow(max.od.matrix),1000)
p.values <- infoScore(max.rate.MI, max.rate.null, nrow(max.rate.matrix),1000)
p.values <- infoScore(carbenicillin.MI, carbenicillin.null, nrow(carbenicillin.matrix),1000)
p.values <- infoScore(tobramycin.MI, tobramycin.null, nrow(tobramycin.matrix),1000)

write

