suppressPackageStartupMessages(library(BGLR))
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(entropy))
suppressPackageStartupMessages(library(gdata))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(TRUE)

parseOGBinary <- function(in.file){
  text.mat <- read.table(in.file, header = T, stringsAsFactors = F)
  num.genomes <- ncol(text.mat)
  num.ogs <- nrow(text.mat)
  genome.list <- colnames(text.mat)
  binary.mat <- matrix(0, nrow = num.ogs, ncol = num.genomes)
  for (i in 1:num.ogs){
    for (j in 1:num.genomes){
      if (text.mat[i,j] != "-"){
        binary.mat[i,j] <- 1
      }
    }
  }
  colnames(binary.mat) <- genome.list
  return(binary.mat)
}
removeHomogeneousFeatures <- function(in.data){
  homogeneous.rows.to.remove <- c()
  for (i in 1:nrow(in.data)){
    if(length(unique(c(in.data[i,]))) == 1) {
      homogeneous.rows.to.remove <- c(homogeneous.rows.to.remove, i)
    }
  }
  n.features <- (length(homogeneous.rows.to.remove))
  if (n.features > 0){
    out.data <- in.data[-homogeneous.rows.to.remove, ]
  } else {
    out.data <- in.data
  }
  ret.list <- list(n.features = n.features, out.data = out.data)
  return(ret.list)
}
mutualInfoEst <- function(matrixData,numGenes) {
  numBins <- 2
  numStrains <- nrow(matrixData)
  mutualMat <- matrix(ncol=numGenes,nrow=numGenes)
  for (i in 1:numGenes) {
    for (j in 1:numGenes) {
      discretVec <- discretize2d(matrixData[,i],matrixData[,j],numBins,numBins)
      mutualMat[i,j] <- suppressWarnings(mi.empirical(discretVec,unit=c("log2")))
    }
  }
  return(mutualMat)
}
nullLagDist <- function(matrixData,numRep,randomizations, stepSize){
  numGenes <- ncol(matrixData)
  distMat <- matrix(ncol=(numGenes^2),nrow=randomizations)
  for (counts in 1:randomizations) {
    randomMatrixData <- matrix(ncol=numGenes,nrow=nrow(matrixData))
    for (i in 1:ncol(matrixData)){
      randomMatrixData[,i] <- sample(matrixData[,i])
    }
    randomMIMat <- lagMIMat(randomMatrixData,numRep, stepSize)
    randomMIVec <- matrix(randomMIMat,ncol=(numGenes^2))
    distMat[counts,] <- randomMIVec
  }
  return(distMat)
}
lagMIMat <- function(matrixData,numRep, stepSize){
  numGenes <- ncol(matrixData)
  numBins <- ceiling(log2(nrow(matrixData)+1))
  step <- numRep * stepSize
  headMatrix <-matrixData[1:(nrow(matrixData)-step),]
  tailMatrix <- matrixData[(step+1):(nrow(matrixData)),]
  lagMatrix <- matrix(ncol=numGenes,nrow=numGenes)
  for (i in 1:numGenes){
    for (j in 1:numGenes){
      lagVec<- suppressWarnings(discretize2d(headMatrix[,i],tailMatrix[,j],numBins,numBins))
      lagMatrix[i,j] <- suppressWarnings(mi.empirical(lagVec,unit=c("log2")))
    }
  }
  return(lagMatrix)
}
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
pangenome.matrix <- 'pa_genomics/genotype/miga-project.ogs'
pheno <- read.table('pa_genomics/phenotype/PA_library_phenotypes_noname.csv', h = T , stringsAsFactors = F, sep = ',')
binary.matrix <- parseOGBinary(in.file = pangenome.matrix)
binary.matrix <- binary.matrix[,pheno$Original.ID]
binary.matrix <- removeHomogeneousFeatures(binary.matrix)$out.data
mutualMat <- mutualInfoEst(matrixData = t(binary.matrix), numGenes = nrow(binary.matrix))

iniMutual <- lagMIMat(matrixData = binary.matrix,numRep = 1,stepSize = 0)
mutualNull <- nullLagDist(matrixData = binary.matrix, ncol(binary.matrix),randomizations = 1000, stepSize = 0)
pValues <- infoScore(iniMutual,mutualNull,numGenes,randomizations)

