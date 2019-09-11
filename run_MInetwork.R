#!/usr/bin/env Rscript

################################################################################
# Name:     run_MInetwork.R
# Author:   Juan C. Castro [jcastro37@gatech.edu]
#           Georgia Institute of Technology
#           
# Version:  1.0
# Date:     07-Aug-2019
# License:  GNU General Public License v3.0.
# ==============================================================================
# 
################################################################################

#======================= 0.0 Install required packages ========================#
personal.lib.path = Sys.getenv("R_LIBS_USER")
if(!file.exists(personal.lib.path))
	            dir.create(personal.lib.path)

packages <- c("entropy", "optparse")
if(any(!(packages %in% installed.packages()))){
	        cat("Please wait a moment! Installing required packages ...\n")
        install.packages(packages[!(packages %in% installed.packages())],
			                          quiet = T, repos="http://cran.rstudio.com/",
						                           lib = personal.lib.path)
	        cat("Required packages installed!\n")
}

#======== 1.0 Load packages, define functions, and initialize variable========#
# 1.1 Load packages ==========================================================#
suppressPackageStartupMessages(library(entropy))
suppressPackageStartupMessages(library(optparse))
# 1.2 Define functions =======================================================#
#' Estimate mutual information (MI) for all pairs of variables in an expression
#' matrix
#'
#' @param matrix.data A matrix with expression data where columns are genes and
#'        rows are time points
#' @param num.genes The number of genes in the matrix
#' @return mutual.mat A matrix with values of mutual information for all pair of
#'        genes 
#' @author Juan C. Castro \email{jcastro37@gatech.edu}
mutualInfoEst <- function(matrix.data,num.genes) {
  mutual.mat <- matrix(ncol=num.genes, nrow=num.genes)
  for (i in 1:num.genes) {
    for (j in 1:num.genes) {
	    if (i>j) {
		    discret.vec <- cbind(matrix.data[,i],matrix.data[,j])
		    mutual.mat[i,j] <- suppressWarnings(mi.empirical(discret.vec,unit=c("log2")))
	    }
    }
  }
  return(mutual.mat)
}
#' Estimate a null distribution of MI for variables in an expression matrix
#'
#' @param matrix.data A matrix with expression data where columns are genes and
#'        rows are time points
#' @param num.genes The number of genes in the matrix
#' @param randomizations Number of times to randomize expression values
#' @return dist.mat A matrix with values of mutual information for a randomized
#'        expression matrix where each column is a linearized version of a random
#'        expression matrix
#' @author Juan C. Castro \email{jcastro37@gatech.edu}
nullInfoDist <- function(matrix.data, num.genes, randomizations){
  ncol.mat <- (num.genes*(num.genes-1))/2 
  dist.mat <- matrix(ncol = ncol.mat, nrow = randomizations)
  for (counts in 1:randomizations) {
    random.matrix.data <- matrix(ncol = ncol(matrix.data), nrow = nrow(matrix.data))
    for (i in 1:ncol(matrix.data)){
      random.matrix.data[,i] <- sample(matrix.data[,i])
    }
    random.mi.mat <- mutualInfoEst(random.matrix.data, num.genes)
    random.mi.vec <- random.mi.mat[lower.tri(random.mi.mat)]
    dist.mat[counts,] <- random.mi.vec
  }
  return(dist.mat)
}
#' Calculate a score for each value of mutual information in an MI matrix
#'
#' @param initial.matrix A matrix with MI values calculated form an expression
#'        matrix
#' @param null.dist.maxtrix A matrix with null values of MI obatined with nullInfoDist
#' @param num.genes The number of genes in the matrix
#' @param randomizations Number of times to randomize expression values
#' @return p.mat A matrix with pValues of for the MI values in initialMatrix
#' @author Juan C. Castro \email{jcastro37@gatech.edu}
infoScore <- function(initial.matrix, null.dist.matrix, num.genes, randomizations){
  mat.size <- (num.genes*(num.genes-1))/2
  sum.dist <- randomizations
  lineal.ini.mat <- initial.matrix[lower.tri(initial.matrix)]
  p.vec <- c()
  for (i in 1:mat.size){
    mi.counts <- sum(null.dist.matrix[,i] >= lineal.ini.mat[i])
    p.score <- mi.counts/sum.dist
    p.vec[i] <- p.score
  }
  p.mat <- matrix(0, ncol = num.genes, nrow = num.genes)
  p.mat[lower.tri(p.mat)] <- p.vec
  return(p.mat)
}
#1.3 Initialize variables ===================================================#
# 1.3.1 Parser variables #
# Get script name
initial.options <- commandArgs(trailingOnly = FALSE)
script.name     <- basename(sub("--file=", "", initial.options[grep("--file=", initial.options)]))
# Process command line arguments
# Create a parser
option_list = list(
  make_option(c("-p", "--pangenome_matrix"), type = "character", default = NULL,
              help ="The OGs matrix that contains the precense abcense values for each gene in each genome.", metavar = "character"),
  make_option(c("-l", "--ogs_list"), type = "character", default = NULL,
              help ="A file with a list with the OG name from which the network is to be built.", metavar = "character"),
  make_option(c("-n", "--num_reps"), type = "numeric", default = 100000,
              help ="Number of repetitions to calculate the null MI distribution.", metavar = "numeric"),
  make_option(c("-o", "--output"), type = "character", default = "network.tsv",
              help ="Name of the output file, the output is stored as a tab separated table of either the OGs matrix (binary or numerical) or the rarefraction curve.", metavar = "character")
  
)
# Add command line arguments
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
opt_names <- names(opt)
if (is.null(opt$pangenome_matrix)){
  print_help(opt_parser)
  err_str <- 'Argument missing "-p" --pangenome_matrix" must be provided.\n'
  stop(err_str, call.=FALSE)
} else if (is.null(opt$ogs_list)){
  print_help(opt_parser)
  err_str <- 'Argument missing "-l" --ogs_list" must be provided.\n'
  stop(err_str, call.=FALSE)
}
# Parse the command line arguments
pangenome.matrix <- opt$pangenome_matrix
ogs.list <- opt$ogs_list
num.reps <- opt$num_reps
output <- opt$output
#===================== 2.0 Format original data as matrix =====================#
# 2.1 Load data ===============================================================#
ogs.df <- read.table(file = pangenome.matrix, header = T, row.names = 1, stringsAsFactors = F)
ogs.list <- read.table(file = ogs.list, header = F, stringsAsFactors = F)
# 2.2 Subset the data to the OGs of the list
sub.df <- ogs.df[ogs.list$V1,]
matrix.data <- t(sub.df)
num.genes <- ncol(matrix.data)
#====================== 3.0 Create the graph based on MI ======================#
# 3.1 Estimate MI and pvalues
initial.mi <- mutualInfoEst(matrix.data = matrix.data, num.genes = num.genes)
null.mi <- nullInfoDist(matrix.data = matrix.data, num.genes = num.genes, randomizations = num.reps)
pvals.mi <- infoScore(initial.matrix = initial.mi, null.dist.matrix = null.mi, num.genes = num.genes, randomizations = num.reps)
# 3.2 Store edges as a data.frame object
edge.list <- data.frame(OG1 = c(), OG2 = c(), MI = c(), pVal = c())
for (i in 1:num.genes){
  for (j in 1:num.genes)
    if(i>j){
      loc.df <- data.frame (OG1 = ogs.list[i,1], OG2 = ogs.list[j,1], MI = initial.mi[i,j], pVal = pvals.mi[i,j])
      edge.list <- rbind(edge.list, loc.df)
    }
}
#===================== 4.0 Write the edge list of a file ======================#
write.table(x = edge.list, file = output, quote = F, sep = '\t', row.names = F)
#==============================================================================#
