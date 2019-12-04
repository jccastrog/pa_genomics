#!/usr/bin/env Rscript

################################################################################
# Name:	    run_BGLR-RFE.R
# Author:   Juan C. Castro [jcastro37@gatech.edu]
#           Georgia Institute of Technology
#           
# Version:  1.15
# Date:     05-Apr-2016
# License:  GNU General Public License v3.0.
# ==============================================================================
# 
################################################################################

#======================= 0.0 Install required packages ========================#
personal.lib.path = Sys.getenv("R_LIBS_USER")
if(!file.exists(personal.lib.path))
  dir.create(personal.lib.path)

packages <- c("optparse", "BGLR", "R.utils")
if(any(!(packages %in% installed.packages()))){
  cat("Please wait a moment! Installing required packages ...\n")
  install.packages(packages[!(packages %in% installed.packages())],
                   quiet = T, repos="http://cran.rstudio.com/",
                   lib = personal.lib.path)
  cat("Required packages installed!\n")
}

#======== 1.0 Load packages, define functions, and initialize variable========#
# 1.1 Load packages ==========================================================#
suppressPackageStartupMessages(library(BGLR))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(R.utils))
# 1.2 Define functions =======================================================#
#' Remove features which have the same value along all genomes
#'
#' @param in.data A matrix with he genotypic data
#' @return ret.list A list with a matrix of genotypic data without redundant features and the number of features removed
#' @author Juan C. Castro \email{jcastro37@gatech.edu}
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
#' Remove features which are present in too few samples (given a percentage cut-off)
#'
#' @param in.data A matrix with the genotypic data
#' @param samples.retain.cutoff Percentage of samples a feature must be present to be included
#' @return ret.list A list with a matrix of genotypic data without sparse features and the number of features removed
#' @author Juan C. Castro \email{jcastro37@gatech.edu}
removeSparseFeatures <- function(in.data, samples.retain.cutoff){
  unique.cols.to.remove <- c()
  n.samples <- ncol(in.data)
  min.n.samples <- floor(n.samples * samples.retain.cutoff)
  for (i in 2:ncol(in.data)) {
    col = in.data[[i]]
    if (length(col[col >= 1]) < min.n.samples){
      unique.cols.to.remove = c(unique.cols.to.remove, i)
    }
  }
  n.features <- length(unique.cols.to.remove)
  if (n.features > 0){
    out.data <- in.data[, -unique.cols.to.remove]
  } else {
    out.data <- in.data
  }
  ret.list <- list(n.features = n.features, out.data = out.data)
  return(ret.list)
}
#' Run BGLR with a set of given parameters
#'
#' @param X A matrix with the genotypic data.
#' @param y A vector with the phenotipic data.
#' @param alg A characters tring describing the model to be used.
#' @param nIter An integer decribing the number iteratons to be used in the MCMC.
#' @param burnIn An interger describing the number of runs tio burn in.
#' @param thin An integer describing the thin parameter of MCNC.
#' @param R2 A numeric value that denotes the expected variance explained by the model. 
#'
#' @return ret.list A list with a matrix of genotypic data without sparse features and the number of features removed
#' @author Juan C. Castro \email{jcastro37@gatech.edu}
runBGLR <- function(X, y, alg, nIter = 1000, burnIn = 100, thin = 3, R2 = 0.5){
  saveAt <- tempdir()
  df0 <- 5
  S0 <- NULL
  weights <- NULL
  if (alg == 'BayesA' | alg == 'BayesB' | alg == 'BayesC' | alg == 'BL' | alg == 'BRR' | alg == 'FIXED'){
    ETA <- list(list(X = X, model = alg))
  } else if (alg == 'RKHS'){
    dimX <- dim(X)
    n <- dimX[1]
    k <- dimX[2]
    X <- as.matrix(X)
    oneM <- matrix(1, nrow = n, ncol = n)
    NX <- X - ((oneM %*% X)*(1/n))
    K <- (t(NX) %*% NX) *(1/n)
    ETA <- list(list(K = K, model = alg))	
  }
  log <- captureOutput({fit_BGLR <- BGLR(y=y,ETA=ETA,nIter=nIter,burnIn=burnIn,thin=thin,saveAt=saveAt,df0=5,S0=S0,weights=weights,R2=R2)})
  #fit_BGLR <- BGLR(y=y,ETA=ETA,nIter=nIter,burnIn=burnIn,thin=thin,saveAt=saveAt,df0=5,S0=S0,weights=weights,R2=R2)
  rem.str <- paste('rm ',saveAt,'*.dat',sep = '')
  system(rem.str)
  return(fit_BGLR)
}
#' Run BGLR with a set of given parameters for a phenotype matrix with more than one phenotypes
#'
#' @param X A matrix with the genotypic data.
#' @param y A vector with the phenotipic data.
#' @param alg A characters tring describing the model to be used.
#' @param nIter An integer decribing the number iteratons to be used in the MCMC.
#' @param burnIn An interger describing the number of runs tio burn in.
#' @param thin An integer describing the thin parameter of MCNC.
#' @param R2 A numeric value that denotes the expected variance explained by the model. 
#'
#' @return ret.list A list with a matrix of genotypic data without sparse features and the number of features removed
#' @author Juan C. Castro \email{jcastro37@gatech.edu}
fitPhenotypes <- function(phenotypeMatrix, X, col.to.remove = 0, alg, nIter, burnIn, thin, R2){
  mse <- data.frame(phenotype = c(), mse = c(), method = c())
  if (col.to.remove > 0){
    phenotypeMatrix <- phenotypeMatrix[,-col.to.remove]
  }
  trait.names <- colnames(phenotypeMatrix)
  n.phenotypes <- dim(phenotypeMatrix)[2]
  phenotype.names <- tail(names(phenotypeMatrix),-1)
  for (i in 2:n.phenotypes){
    y <- phenotypeMatrix[,i]
    fit.sm  <- runBGLR(X, y,  alg, nIter, burnIn, thin, R2)
    local.mse <- data.frame(phenotype = phenotype.names[i-1], mse = mean((fit.sm$y - fit.sm$yHat)^2), method = alg)
    mse <- rbind(mse, local.mse)
    cat(paste(i-1,' of ',n.phenotypes-1,' phenotypes processed...\n',sep=''))
  }
  return(mse)
}
#' Select random balanced folds
#'
#' @param X  A matrix with the genotypic data.
#' @param num.folds The number of folds to create in the data
#'
#' @return partition.matrix A matrix with the indexes for the randoms folds in the X 
#'         matrix (The same indexes are used to select the folds of y) with num.folds.
#'         number of rows.
#' @author Juan C. Castro \email{jcastro37@gatech.edu}
selectFolds <- function(X, num.folds = 10){
  num.genomes <- dim(X)[1]
  fold.size <- num.genomes/num.folds
  big.fold.size <- ceiling(fold.size)
  small.fold.size <- floor(fold.size)
  num.bigs <- num.genomes%%num.folds
  rand.perm <- sample(seq(1,num.genomes), replace = F)
  partition.matrix <- matrix(0, ncol = big.fold.size, nrow = num.folds)
  k = 1
  while (k <= num.bigs){
    begin <- (k * big.fold.size) - (big.fold.size - 1)
    end <- k * big.fold.size
    partition.matrix[k,] <- rand.perm[begin:end]
    k <- k + 1
  }
  while (k <= num.folds){
    begin <- (k * small.fold.size) - (small.fold.size - 1)
    end <- k * small.fold.size
    partition.matrix[k,-big.fold.size] <- rand.perm[begin:end]
    k <- k + 1
  }
  return(partition.matrix)
}
#' Execute BGLR-RFE
#'
#' @param X  A matrix with the genotypic data.
#' @param phenotype.data A matrix with phenotypic data 
#' @param folds.matrix A matrix with the selected cross validatiion folds
#' @param trait The column in the phenotype.data matrix that is to be predicted
#' @param alg A characters tring describing the model to be used.
#' @param n.iter An integer decribing the number iteratons to be used in the MCMC.
#' @param burn.in An interger describing the number of runs tio burn in.
#' @param thin An integer describing the thin parameter of MCNC.
#' @param R2 A numeric value that denotes the expected variance explained by the model. 
#'
#' @return partition.matrix A matrix with the indexes for the randoms folds in the X 
#'         matrix (The same indexes are used to select the folds of y) with num.folds.
#'         number of rows.
#' @author Juan C. Castro \email{jcastro37@gatech.edu}
runBGLR_RFE <- function(folds.matrix, X, phenotype.data, trait, alg, n.iter, burn.in, thin, R2){
  n.features <- ncol(X)
  num.folds <- nrow(folds.matrix)
  rank.matrix <- data.frame(rank=seq(n.features,1))
  mse.table <- data.frame(fold = c(), test = c(), all = c())
  og.names <- colnames(X)
  params.effects <- data.frame(row.names = og.names)
  params.sd <- data.frame(row.names = og.names)
  for (i in 1:num.folds){
    k <- 1
    X.fold <- X
    n.features <- dim(X.fold)[2]
    ranks <- data.frame(og = c())
    cat(paste('Processing fold ',i,' of ',num.folds,'...\n',sep = ''))
    pheno.fold <- phenotype.data
    pheno.fold[folds.matrix[i,],] <- NA
    y.fold <- as.numeric(pheno.fold[[trait]])
    y.test <- as.numeric(phenotype.data[[trait]])
    while (k <= n.features) {
      if (k == 1){
        fit.BGLR <- runBGLR(X.fold, y.fold, alg = algorithm)
        mse.fold <- mean((fit.BGLR$y - fit.BGLR$yHat)^2, na.rm = TRUE)
        mse.test <- mean((y.test[folds.matrix[i,]] - predict(fit.BGLR)[folds.matrix[i,]])^2)
        mse.all <- mean((y.test - predict(fit.BGLR))^2)
        param.fold <- fit.BGLR$ETA[[1]]$b
        param.sd <- fit.BGLR$ETA[[1]]$SD.b
        params.effects[, i] <- cbind(params.effects, param.fold[match(rownames(params.effects),names(param.fold))])
        params.sd[, i] <- cbind(params.sd, param.sd[match(rownames(params.sd),names(param.sd))])
        sorted.param.fold <- sort(abs(param.fold), decreasing = F)
        names.param.fold <- names(sorted.param.fold)
        k <- k + 1
        temp.mse <- data.frame(fold = mse.fold, test = mse.test, all = mse.all)
        mse.table <- rbind(mse.table,temp.mse)
        #Remove last feature
        rem.feat <- as.character(tail(names.param.fold, n = 1))
        rem.index <- match(rem.feat, colnames(X.fold))
        temp.rank <- data.frame(og = rem.feat)
        ranks <- rbind(ranks, temp.rank)
        X.fold <- X.fold[,-rem.index]
      } else {
        fit.BGLR <- runBGLR(X.fold, y.fold, alg = algorithm, nIter = n.iter, burnIn = burn.in, thin = thin, R2 = R2)
        param.fold <- fit.BGLR$ETA[[1]]$b
        sorted.param.fold <- sort(param.fold, decreasing = F)
        names.param.fold <- names(sorted.param.fold)
        rem.feat <- as.character(tail(names.param.fold, n = 1))
        rem.index <- match(rem.feat, colnames(X.fold))
        if (dim(X.fold)[2] > 2){
          temp.rank <- data.frame(og = rem.feat)
          ranks <- rbind(ranks, temp.rank)
          X.fold <- X.fold[,-rem.index]
        } else if(dim(X.fold)[2] == 2){
          keep.feat <- as.character(head(names.param.fold, n = 1))
          temp.rank <- data.frame(og = c(rem.feat, keep.feat))
          ranks <- rbind(ranks, temp.rank)
          k <- n.features + 1
        }
        k <- k + 1
      }
      #cat(paste('Feature ', k, ' of ', n.features,'\n'))
    }
    rank.matrix <- cbind(rank.matrix,ranks)
    rank.matrix <- rank.matrix[order(rank.matrix$rank),]
    rank.matrix <- averageRankMatrix(rank.matrix)
    mean.effects <- apply(params.effects, 1, mean)
    mean.sd <- apply(params.sd, 1, mean)
    mean.t.ratio <- mean.effects/mean.sd
    mean.p.value <- pt(-abs(mean.t.ratio), df = 5)
    adjust.p.value <- p.adjust(mean.p.value, method = "fdr", n = length(mean.p.value))
    reg.effects <- cbind(mean.effects, mean.sd, mean.p.value, adjust.p.value)
    ret.list <- list(mse = mse.table, rank.matrix = rank.matrix, effects = reg.effects)
  }
  return(ret.list)
}
#' Average ranks
#'
#' @param ranked.matrix  A matrix a coulnm for the rank of each OG in each fold
#'
#' @return sorted.rank A table with the rank for each OG 
#' @author Juan C. Castro \email{jcastro37@gatech.edu}
averageRankMatrix <- function(ranked.matrix){
  num.features <- nrow(ranked.matrix)
  num.folds <- ncol(ranked.matrix) - 1
  total.ranks <- matrix(NA, nrow = num.features, ncol = num.folds)
  rownames(total.ranks) <- unique(ranked.matrix[,2])
  for (i in 1:num.folds){
    temp.lst <- as.character(ranked.matrix[,i + 1])
    for (j in 1:num.features){
      og <- temp.lst[j]
      total.ranks[og,i] <- j
    }
  }
  mean.ranks <- apply(total.ranks, 1, mean)
  sorted.ogs <- names(sort(mean.ranks))
  final.ranks <- data.frame(rank = seq(1, num.features), og = sorted.ogs)
  return(final.ranks)
}
# 1.3 Initialize variables ===================================================#
# 1.3.1 Parser variables #
# Get script name
initial.options <- commandArgs(trailingOnly = FALSE)
script.name     <- basename(sub("--file=", "", initial.options[grep("--file=", initial.options)]))
# Process command line arguments
# Create a parser
option_list = list(
  make_option(c("-g", "--genotype_matrix"), type = "character", default = NULL,
              help ="The genotyoe matrix file.", metavar = "character"),
  make_option(c("-p", "--phenotype_matrix"), type = "character", default = "NULL",
              help ="A file containing the phenotype values, each column is a different phenotype.", metavar = "character"),
  make_option(c("-q", "--quirk"), type = "integer", default = "NULL",
              help ="The column in the phenotype file that contains the trait to be predicted.", metavar = "integer"),
  make_option(c("-a", "--algorithm"), type = "character", default = "BayesA",
              help ="The algorithms to be used (BayesA, BayesB, BayesC, RKHS, BL). [default= %default]", metavar = "character"),
  make_option(c("-f", "--num_folds"), type = "integer", default = "6",
              help ="Number of folds for the corss validation procedure [default= %default].", metavar = "integer"),
  make_option(c("-i", "--num_iter"), type = "integer", default = 1000,
              help ="Number of iterations of the MCMC. [default= %default]", metavar = "integer"),
  make_option(c("-b", "--burn_in"), type = "integer", default = 100,
              help ="Burn-in value for the MCMC. [default= %default]", metavar = "integer"),
  make_option(c("-t", "--thin"), type = "integer", default = 1,
              help ="Thining value of the MCMC. [default= %default]", metavar = "integer"),
  make_option(c("-r", "--r2"), type = "numeric", default = 0.5,
              help ="The expected variance to be explained by the model. [default= %default]", metavar = "numeric"),
  make_option(c("-s", "--samp_perc"), type = "numeric", default = 0.01,
              help ="Percentage of samples in which a feature must be present to be retained in the set. [default= %default]", metavar = "numeric"),
  make_option(c("-c", "--remove_col"), type = "integer", default = 0,
              help = "Column number to be removed if any  [default= %default]", metavar = "integer"),
  make_option(c("-n", "--names_col"), type = "integer", default = 0,
              help = "In the phenotype matrix the column number with the genome with the IDs", metavar = "integer"),
  make_option(c("-o", "--output"), type = "character", default = "rank_list.txt",
              help ="Output filename [default= %default]", metavar = "character")
);
# Add command line arguments
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
opt_names <- names(opt)
if (is.null(opt$genotype_matrix)){
  print_help(opt_parser)
  err_str <- 'Argument missing "--genotype_matrix" must be provided.\n'
  stop(err_str, call.=FALSE)
}
# Parse the command line arguments
genotype.matrix <- opt$genotype_matrix
phenotype.matrix <- opt$phenotype_matrix
quirk <- opt$quirk
algorithm <- opt$algorithm
num.folds <- opt$num_folds
num.iter <- opt$num_iter
burn.in <- opt$burn_in
thin <- opt$thin
R2 <- opt$r2
samp.perc.retain <- opt$samp_perc
names.col <- opt$names_col
remove.col <- opt$remove_col
output <- opt$output
# 1.3.2 Load data #
cat('Loading data... ')
genotype.data <- read.table(genotype.matrix, h = T, sep = '\t', stringsAsFactors = F, row.names = 1)
phenotype.data <- read.table(phenotype.matrix, h = T, sep = ',', stringsAsFactors=F)
names.phenotype <- phenotype.data[,names.col]
if (!any(is.na(remove.col))) {
  phenotype.data <- phenotype.data[-remove.col]
}
n.genomes <- ncol(genotype.data)
genotype.data <- genotype.data[names.phenotype]
phenotype.data <- phenotype.data[-names.col]
rownames(phenotype.data) <- names.phenotype
data.dim <- dim(genotype.data)
n.phenotypes <- nrow(phenotype.data)
n.traits <- ncol(phenotype.data)
cat('Done!\n')
cat(paste(n.genomes,' genomes and ',n.phenotypes,' phenotypes loaded!\n'))
#============== 2.0 Clean data and prepare for cross-validation ==============#
# 2.1 Remove redundant and sparse data =======================================# 
cat('Removing features with same values across all samples... ')
non.homogeneous.data <- removeHomogeneousFeatures(genotype.data)
n.features.removed <- non.homogeneous.data$n.features
non.sparse.data <- removeSparseFeatures(t(non.homogeneous.data$out.data), samp.perc.retain)
clean.data <- non.sparse.data$out.data
n.features <- ncol(clean.data)
n.features.removed <- n.features.removed + non.sparse.data$n.features
n.ogs <- dim(clean.data)[2]
cat('Done!\n')
flush.console()
cat(paste(n.features.removed, " features removed\n"), sep = "")
flush.console()
# 2.2 Select folds ===========================================================#
cat('Selecting balanced folds for cross-validation...')
folds.matrix <- selectFolds(clean.data, num.folds)
cat('Done!\n')
flush.console()
#====================== 3.0 Run BGLR-RFE for each fold ======================#
# 3.1 Execute BGLR-RFE ======================================================#
cat('Executing BGLR-RFE ...\n')
rfe.list <- runBGLR_RFE(folds.matrix = folds.matrix, X= clean.data, phenotype.data = phenotype.data, trait = quirk-2, alg = algorithm, n.iter = num.iter, burn.in = burn.in, thin = thin, R2 = R2)
# 3.2 Save the results in output files ======================================#
write.table(x = rfe.list$mse, file = paste(output,"_mse.tsv", sep = ''), sep = '\t', quote = F, row.names = T)
write.table(x = rfe.list$rank.matrix, file = paste(output,"_ranks.tsv", sep = ''), sep = '\t', quote = F, row.names = T)
write.table(x = rfe.list$effects, file = paste(output,"_effects.tsv", sep = ''), sep = '\t', quote = F, row.names = T)
cat(paste('Done, results can be found in ',output,'\n',sep=''))
#=============================================================================#
#rfe <- runBGLR_RFE(folds.matrix = folds.matrix, X= clean.data[,1:1000], phenotype.data = phenotype.data, trait = 1, alg = algorithm, n.iter = num.iter, burn.in = burn.in, thin = thin, R2 = 0.5)
#rfe.biofilm.od <- runBGLR_RFE(folds.matrix = folds.matrix, X= clean.data, phenotype.data = phenotype.data, trait = 1, alg = algorithm, n.iter = num.iter, burn.in = burn.in, thin = thin, R2 = 0.5)
#rfe.biofilm.eps <- runBGLR_RFE(folds.matrix = folds.matrix, X= clean.data, phenotype.data = phenotype.data, trait = 2, alg = algorithm, n.iter = num.iter, burn.in = burn.in, thin = thin, R2 = 0.5)
#rfe.fraction.biofilm <- runBGLR_RFE(folds.matrix = folds.matrix, X= clean.data, phenotype.data = phenotype.data, trait = 3, alg = algorithm, n.iter = num.iter, burn.in = burn.in, thin = thin, R2 = 0.5)
#rfe.planktonic.od <- runBGLR_RFE(folds.matrix = folds.matrix, X= clean.data, phenotype.data = phenotype.data, trait = 4, alg = algorithm, n.iter = num.iter, burn.in = burn.in, thin = thin, R2 = 0.5)
#rfe.planktonic.od.variance <- runBGLR_RFE(folds.matrix = folds.matrix, X= clean.data, phenotype.data = phenotype.data, trait = 5, alg = algorithm, n.iter = num.iter, burn.in = burn.in, thin = thin, R2 = 0.5)
#rfe.max.od <- runBGLR_RFE(folds.matrix = folds.matrix, X= clean.data, phenotype.data = phenotype.data, trait = 6, alg = algorithm, n.iter = num.iter, burn.in = burn.in, thin = thin, R2 = 0.5)
#rfe.max.rate <- runBGLR_RFE(folds.matrix = folds.matrix, X= clean.data, phenotype.data = phenotype.data, trait = 7, alg = algorithm, n.iter = num.iter, burn.in = burn.in, thin = thin, R2 = 0.5)
#rfe.elastase <- runBGLR_RFE(folds.matrix = folds.matrix, X= clean.data, phenotype.data = phenotype.data, trait = 8, alg = algorithm, n.iter = num.iter, burn.in = burn.in, thin = thin, R2 = 0.5)
#rfe.pyochelin <- runBGLR_RFE(folds.matrix = folds.matrix, X= clean.data, phenotype.data = phenotype.data, trait = 9, alg = algorithm, n.iter = num.iter, burn.in = burn.in, thin = thin, R2 = 0.5)
#rfe.pyoverdin <- runBGLR_RFE(folds.matrix = folds.matrix, X= clean.data, phenotype.data = phenotype.data, trait = 10, alg = algorithm, n.iter = num.iter, burn.in = burn.in, thin = thin, R2 = 0.5)
#rfe.sa.virulence <- runBGLR_RFE(folds.matrix = folds.matrix, X= clean.data, phenotype.data = phenotype.data, trait = 11, alg = algorithm, n.iter = num.iter, burn.in = burn.in, thin = thin, R2 = 0.5) 
#rfe.control.auc <- runBGLR_RFE(folds.matrix = folds.matrix, X= clean.data, phenotype.data = phenotype.data, trait = 12, alg = algorithm, n.iter = num.iter, burn.in = burn.in, thin = thin, R2 = 0.5)
#rfe.carbenicillin.auc <- runBGLR_RFE(folds.matrix = folds.matrix, X= clean.data, phenotype.data = phenotype.data, trait = 13, alg = algorithm, n.iter = num.iter, burn.in = burn.in, thin = thin, R2 = 0.5)
#rfe.tobramycin.auc <- runBGLR_RFE(folds.matrix = folds.matrix, X= clean.data, phenotype.data = phenotype.data, trait = 14, alg = algorithm, n.iter = num.iter, burn.in = burn.in, thin = thin, R2 = 0.5)



