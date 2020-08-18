#!/usr/bin/env Rscript

################################################################################
# Name:	    mutualInfoHeuristic.R
# Author:   Juan C. Castro [jcastro37@gatech.edu]
#           Georgia Institute of Technology
#           
# Version:  1.0
# Date:     18-Aug-2020
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
# 1.2 Define functions =======================================================#
#' Estimate mutual information (MI) for all pairs of variables in an expression
#' matrix
#'
#' @param matrixData A matrix with binary values where the colums are genomes
#'        and the rows are orthologous groups
#' @return mutualMat A matrix with values of mutual information for all pair of
#'        genes 
#' @author Juan C. Castro \email{jcastro37@gatech.edu}
mutualInfoEst <- function(matrixData) {
  numGenes <- nrow(matrixData)
  numBins <- 2
  mutualMat <- matrix(ncol=numGenes,nrow=numGenes)
  for (i in 1:numGenes) {
    for (j in 1:numGenes) {
      if (i>j){
        numVarsI <- length(table(matrixData[i,]))
        numVarsJ <- length(table(matrixData[j,]))
        if (numVarsI >1 & numVarsJ >1){
          discretVec <- discretize2d(matrixData[i,],matrixData[j,],numBins,numBins)
          mutualMat[i,j] <- suppressWarnings(mi.empirical(discretVec,unit=c("log2")))
        } else {
          mutualMat[i,j] <- 0
        }
      }
    }
  }
  return(mutualMat)
}
#' Estimate a null distribution of MI for variables in an expression matrix
#'
#' @param numGenomes The number of genomes in the matrix for which the heuristic
#'        is being calculated
#' @param randomizations Number of times to randomize expression values
#' @return MIvals A vector with a sample of posible values of mutual information
#'        for a given matrix size
#' @author Juan C. Castro \email{jcastro37@gatech.edu}
nullInfoDist <- function(numGenomes,randomizations){
  MIvals <- c()
  for (counts in 1:randomizations){
    vecI <- round(runif(n = numGenomes, min = 0, max = 1))
    vecJ <- round(runif(n = numGenomes, min = 0, max = 1))
    numVarsI <- length(table(vecI))
    numVarsJ <- length(table(vecJ))
    if (numVarsI >1 & numVarsJ >1){
      discretVec <- discretize2d(vecI, vecJ, 2, 2)
      locMI <- suppressWarnings(mi.empirical(discretVec,unit=c("log2")))
    } else {
      locMI <- 0
    }
    MIvals <- c(MIvals, locMI)
  }
  return(MIvals)
}
#' Calculate a score for each value of mutual information in an MI matrix
#'
#' @param initialMatrix A matrix with MI values calculated form an expression
#'        matrix
#' @param nullDistVec A vector with a sample of posible values of mutual information
#'        for a given matrix size
#' @param numGenes The number of genes in the matrix
#' @return pMat A matrix with pValues of for the MI values in initialMatrix
#' @author Juan C. Castro \email{jcastro37@gatech.edu}
infoScore <- function(initialMatrix,nullDistVec){
  numOGS <-  ncol(initialMatrix)
  sumDist <- length(nullDistVec)
  pMat <- matrix(0, ncol = numOGS , nrow = numOGS)
  for (i in 1:numOGS){
    for (j in 1:numOGS){
      if (i>j){
        pScore <- sum(nullDistVec >= initialMatrix[i,j])/sumDist
        pMat[i,j] <- pScore
      }
    }
  }
  return(pMat)
}
#1.3 Initialize variables ===================================================#
# 1.3.1 Parser variables #
# Get script name
initial.options <- commandArgs(trailingOnly = FALSE)
script.name     <- basename(sub("--file=", "", initial.options[grep("--file=", initial.options)]))
# Process command line arguments
# Create a parser
option_list = list(
  make_option(c("-g", "--pangenome_matrix"), type = "character", default = NULL,
              help ="The OGs matrix in which each row is an OG and each column is a strain.", metavar = "character"),
  make_option(c("-l", "--list_ogs"), type = "logical", default = FALSE,
              help ="A lilst with the OGs selected as relevant features.", metavar = "logical"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help ="Name of the output file, the output is stored as a tab separated table of either the OGs matrix (binary or numerical) or the rarefraction curve.", metavar = "character")
  
)
# Add command line arguments
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
opt_names <- names(opt)
if (is.null(opt$pangenome_matrix)){
  print_help(opt_parser)
  err_str <- 'Argument missing "-g" --pangenome_matrix" must be provided.\n'
  stop(err_str, call.=FALSE)
} else if (is.null(opt$output)){
  print_help(opt_parser)
  err_str <- 'Argument missing "-o" --output" must be provided.\n'
  stop(err_str, call.=FALSE)
}
# Parse the command line arguments
pangenome.matrix <- opt$pangenome_matrix
if (exists("opt$list_ogs")){
  list.ogs <- opt$list_ogs
}
output <- opt$output
#Global variables
big.edge.list <- data.frame(g1 = c(), g2 = c(), pval = c()) 
#============== 2.0 Load and process the data ==============#
# 2.1 Load the data ========================================#
pangenome.df <- read.table(pangenome.matrix, h = T ,
                           stringsAsFactors = F, sep = '\t', row.names = 1)
if (exists("list.ogs")){
ogs <- read.table(list.ogs, h = F,
                  stringsAsFactors = F)
}
# 2.2 Process the data =====================================#
if (exists("ogs")){
  sub.pangenome <- subset(pangenome.df, rownames(pangenome.df)%in%ogs)
} else {
  sub.pangenome <- pangenome.df
}
#============== 3.0 Calculate MI values and p.values respectively ==============#
initial.MI <- mutualInfoEst(matrixData = sub.pangenome)
null.MI <- nullInfoDist(numGenomes = ncol(sub.pangenome), randomizations = 1000000)
p.values < - infoScore(initialMatrix = initial.MI, nullDistVec = null.MI)
#============== 4.0 Parse the p.values as edge list ==============#
i.names <- row.names(p.values)
j.names <- col.names(p.values)
for (i in 1:nrow(p.values)){
  for (j in 1:ncol(p.values)){
    loc.edge.list <- data.frame(g1 = i.names[i], g2 = j.names[j], pval = p.values[i,j])     
    big.edge.list <- rbind(big.edge.list, loc.edge.list)
  }
}
write.table(x = big.edge.list, file = output, quote = F, sep = '\t')
#==========================================================================================#
small.edge.list <- subset(big.edge.list, big.edge.list$pval < 0.01)
local.net <- graph_from_data_frame(d = small.edge.list[,1:2], directed = F)
#==========================================================================================#