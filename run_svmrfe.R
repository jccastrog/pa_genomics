#!/usr/bin/env Rscript

################################################################################
# Name:	    run_svmrfe.R
# Author:   Juan C. Castro [jcastro37@gatech.edu]
#           Georgia Institute of Technology
#           
# Version:  1.0
# Date:     23-Nov-2018
# License:  GNU General Public License v3.0.
# ==============================================================================
# 
################################################################################


#======================= 0.0 Install required packages ========================#
personal.lib.path = Sys.getenv("R_LIBS_USER")
if(!file.exists(personal.lib.path))
  dir.create(personal.lib.path)

packages <- c("e1071", "optparse")
if(any(!(packages %in% installed.packages()))){
  cat("Please wait a moment! Installing required packages ...\n")
  install.packages(packages[!(packages %in% installed.packages())],
                   quiet = T, repos="http://cran.rstudio.com/",
                   lib = personal.lib.path)
  cat("Required packages installed!\n")
}
#======== 1.0 Load packages, define functions, and initialize variable========#
# 1.1 Load packages ==========================================================#
suppressPackageStartupMessages(library(optparse, quietly = TRUE))
suppressPackageStartupMessages(library(e1071, quietly = TRUE))
# 1.2 Define functions =======================================================#
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
  n.samples <- nrow(in.data)
  min.n.samples <- floor(n.samples * samples.retain.cutoff)
  for (i in 1:ncol(in.data)) {
    col = in.data[,i]
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
    if (num.bigs != 0){
      partition.matrix[k,-big.fold.size] <- rand.perm[begin:end]
    } else {
      partition.matrix[k,] <- rand.perm[begin:end]
    }
    k <- k + 1
  }
  return(partition.matrix)
}
svmRFE <- function(X, cls, folds.matrix, k = 1){
  n.features <- ncol(X)
  num.folds <- nrow(folds.matrix)
  rank.matrix <- data.frame(rank=seq(n.features,1))
  og.names <- colnames(X)
  for (i in 1:num.folds){
    k <- 1
    X.fold <- X
    ranks <- data.frame(og = c())
    cat(paste('Processing fold ',i,' of ',num.folds,'...\n',sep = ''))
    y.fold <- cls
    #y.fold[folds.matrix[i,]] <- NA
    while (k <= n.features) {
      if (k == 1){
        weights <- as.data.frame(getWeights(y.fold, X.fold))
        sorted.weights <- sort(abs(weights), decreasing = F)
        names.weights <- names(sorted.weights)
        k <- k + 1
        rem.feat <- as.character(head(names.weights, n = 1))
        rem.index <- match(rem.feat, colnames(X.fold))
        temp.rank <- data.frame(og = rem.feat)
        ranks <- rbind(ranks, temp.rank)
        X.fold <- X.fold[,-rem.index]
      } else {
        weights <- as.data.frame(getWeights(y.fold, X.fold))
        sorted.weights <- sort(abs(weights), decreasing = F)
        names.weights <- names(sorted.weights)
        rem.feat <- as.character(tail(names.weights, n = 1))
        rem.index <- match(rem.feat, colnames(X.fold))
        if (dim(X.fold)[2] > 2){
          temp.rank <- data.frame(og = rem.feat)
          ranks <- rbind(ranks, temp.rank)
          X.fold <- X.fold[,-rem.index]
        } else if(dim(X.fold)[2] == 2){
          keep.feat <- as.character(head(names.weights, n = 1))
          temp.rank <- data.frame(og = c(rem.feat, keep.feat))
          ranks <- rbind(ranks, temp.rank)
          rem.feat <- as.character(head(names.weights, n = 1))
          rem.index <- match(rem.feat, colnames(X.fold))
          k <- n.features + 1
        }
        k <- k + 1
      }
    }
    rank.matrix <- cbind(rank.matrix,ranks)
  }
  rank.matrix <- rank.matrix[order(rank.matrix$rank),]
  rank.matrix <- averageRankMatrix(rank.matrix)
  return <- rank.matrix
}

getWeights <- function(y, X) {
  # Fit a linear SVM model and obtain feature weights
  train.data = y
  svmModel = svm(X, train.data, 
                 cost=100, cachesize=1000,
                 scale=F, type="C-classification", 
                 kernel="linear")
  
  return(t(svmModel$coefs) %*% svmModel$SV)
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
  make_option(c("-i", "--num_iter"), type = "integer", default = 1000,
              help ="Number of iterations of the MCMC. [default= %default]", metavar = "integer"),
  make_option(c("-s", "--samp_perc"), type = "numeric", default = 0.01,
              help ="Percentage of samples in which a feature must be present to be retained in the set. [default= %default]", metavar = "numeric"),
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
num.folds <- opt$num_folds
samp.perc.retain <- opt$samp_perc
names.col <- opt$names_col
output <- opt$output
###================FOR DEBUGGING PURPOUSES ONLY===============###
#################################################################
genotype.matrix <- 'F:/JuanCamilo/GaTech/BrownLab/projects/pa_genomics/genotype/pa_genomes.ogs.numerical.codes'
#genotype.matrix <- '/media/jccastrog/Data/JuanCamilo/GaTech/BrownLab/projects/pa_genomics/genotype/pa_genomes.ogs.numerical.codes'
phenotype.matrix <- 'F:/JuanCamilo/GaTech/BrownLab/projects/pa_genomics/phenotype/PA_library_phenotypes_noname.csv'
#phenotype.matrix <- '/media/jccastrog/Data/JuanCamilo/GaTech/BrownLab/projects/pa_genomics/phenotype/PA_library_phenotypes_noname.csv'
num.folds <- 10
samp.perc.retain <- 0.05
names.col <- 1
# 1.3.2 Load data #
cat('Loading data... ')
genotype.data <- read.table(genotype.matrix, h = T, sep = '\t', stringsAsFactors = F, row.names = 1)
phenotype.data <- read.table(phenotype.matrix, h = T, sep = ',', stringsAsFactors=F)
names.phenotype <- phenotype.data[,names.col]
n.genomes <- ncol(genotype.data)
genotype.data <- genotype.data[names.phenotype]
phenotype.data <- phenotype.data[-names.col]
rownames(phenotype.data) <- names.phenotype
data.dim <- dim(genotype.data)
n.phenotypes <- nrow(phenotype.data)
n.traits <- ncol(phenotype.data)
cat('Done!\n')
phenotype.data$classes <- NA
for (i in 1:length(phenotype.data$classes)){
  if (phenotype.data$env.source[i] == "cf"){
    phenotype.data$classes[i] <- 1
  } else {
    phenotype.data$classes[i] <- 0
  }
}
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
folds.matrix <- selectFolds(clean.data, 1)
cat('Done!\n')
flush.console()
rank.matrix <- svmRFE(X = clean.data, cls = phenotype.data$classes, folds.matrix = folds.matrix )

top.rank <- rank.matrix[1:10,]
sub.mat <- as.matrix(data.frame(OG_57 = clean.data[,"OG_57"], OG_72 = clean.data[,"OG_72"], OG_1328 = clean.data[,"OG_1328"],
                      OG_1354 = clean.data[,"OG_1354"], OG_1699 = clean.data[,"OG_1699"], OG_1899 = clean.data[,"OG_1889"],
                      OG_2527 = clean.data[,"OG_2527"], OG_2559 = clean.data[,"OG_2559"], OG_2745 = clean.data[,"OG_2745"],
                      OG_2746 = clean.data[,"OG_2746"]))



sub.mat <- clean.data[,top.rank$og]
my.palette <- colorRampPalette(colors = c("orange", "#00ffff", "blue"))(n=3)
pdf('plots/heatmap_SVM_CF.pdf', width = 7, height = 7)
tiff('plots/heatmap_SVM_CF.tiff',res = 300, width = 3000, height = 3000)
heatmap.2(sub.mat, trace = "none", dendrogram = "row", density.info = "none", col = my.palette)
graphics.off()


#SHEYDA-SVM======================================================================================
removeNA <- function(X, samp.perc.retain){
  samp.tot <- ncol(X)
  remove.threshold <- (1-samp.perc.retain) * samp.tot
  rows.to.remove <- c()
  for (i in 1:nrow(X)){
    loc.row <- X[i,]
    row.na <- loc.row[is.na(loc.row)]
    num.na <- length(row.na)
    if (num.na >= remove.threshold){
      rows.to.remove <- c(rows.to.remove, i)
    }
  }
  new.tab <- X[-rows.to.remove,]
  return(new.tab)
}

X <- t(non.sparse.data$out.data)
clean.data <- removeNA(t(non.sparse.data$out.data), samp.perc.retain = 0.05)
rank.matrix <- svmRFE(X = t(clean.data), cls = phenotype.data$Piperacillin_Tazobactam, folds.matrix = folds.matrix )
rank.matrix <- svmRFE(X = t(clean.data), cls = phenotype.data$Ceftazidim, folds.matrix = folds.matrix )

svmModel = svm(x =t(clean.data), y = c(phenotype.data$Piperacillin_Tazobactam), cost=100, cachesize=1000, scale=F, type="C-classification", kernel="linear")
weights <- t(svmModel$coefs) %*% matrix(as.numeric(svmModel$SV), ncol = 13)
plot(c(1:13), weights)

svmModel = svm(x =t(clean.data), y = c(phenotype.data$Ceftazidim), cost=100, cachesize=1000, scale=F, type="C-classification", kernel="linear")
weights <- t(svmModel$coefs) %*% matrix(as.numeric(svmModel$SV), ncol = 13)
plot(c(1:13), weights)

markers <- c(556486, 4173225)

heatmap.2clean.data[c(3,7),]