#!/usr/bin/env Rscript

################################################################################
# Name:	    getGeneCorRanks.R
# Author:   Juan C. Castro [jcastro37@gatech.edu]
#           Georgia Institute of Technology
#           
# Version:  1.0
# Date:     27-Dec-2018
# License:  GNU General Public License v3.0.
# ==============================================================================
# 
################################################################################

#======================= 0.0 Install required packages ========================#
personal.lib.path = Sys.getenv("R_LIBS_USER")
if(!file.exists(personal.lib.path))
  dir.create(personal.lib.path)

packages <- c("ape", "optparse", "phylosignal")
if(any(!(packages %in% installed.packages()))){
  cat("Please wait a moment! Installing required packages ...\n")
  install.packages(packages[!(packages %in% installed.packages())],
                   quiet = T, repos="http://cran.rstudio.com/",
                   lib = personal.lib.path)
  cat("Required packages installed!\n")
}
#======== 1.0 Load packages, define functions, and initialize variable========#
# 1.1 Load packages ==========================================================#
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(phylosignal))
# 1.2 Define functions =======================================================#
#' Take a phenotype table and a tree and make a
#'
#' @param tree.file A file in nwk format with an estimated phylogeny of an OG
#' @param phenotype.table A table where the row names are the strain names and
#'                        the columns are phenotypic traits, each entry is a
#'                        measurement for a phenotype in a strain
#' @return traits.signal A lilst with stats and p values for the Cmean value
#'                       for each trait in phenotype.table
#' @author Juan C. Castro \email{jcastro37@gatech.edu}
getPhyloCmean <- function(tree.file, phenotype.table){
  og.tree <- read.tree(tree.file)
  sub.phenotype <- phenotype.table[og.tree$tip.label,]
  temp.pheno <- tempfile()
  write.table(x = sub.phenotype, file = temp.pheno, sep = '\t', row.names = T, quote = F)
  local.p4d <- read.p4d(phylo.file = tree.file, data.file = temp.pheno)
  traits.signal <- phyloSignal(p4d = local.p4d, methods = "Cmean")
  return(traits.signal)
}
#' Get a directory of tree files and get the matrices for the statistic 
#' Cmean and P values. 
#'
#' @param in.file Path to the file containing the matrix 
#' @return og.matrix A matrix with the information regarding the OGs in binary
#'                   terms where 1 is precense of the OG and 0 mean abscence 
#' @author Juan C. Castro \email{jcastro37@gatech.edu}
getMatrices <- function(tree.dir, phenotype.table, tree.extension, min.num.strains){
  tree.files <- dir(tree.dir, tree.extension)
  num.ogs <- length(tree.files)
  stats.matrix <- matrix(0, nrow = num.ogs, ncol = ncol(phenotype.table))
  pval.matrix <- matrix(0, nrow = num.ogs, ncol = ncol(phenotype.table))
  colnames(stats.matrix) <- colnames(phenotype.table)
  colnames(pval.matrix) <- colnames(phenotype.table)
  og.names <- c()
  k.og <-1
  for (i in tree.files){
    loc.tree.file <- paste(tree.dir,"/", i, sep = "")
    og.loc.name <- gsub("RAxML_bestTree.", "", i)
    og.loc.name <- gsub(".out", "", og.loc.name)
    og.names <- c(og.names, og.loc.name)
    og.tree <- read.tree(loc.tree.file)
    if (length(og.tree$tip.label) >= min.num.strains){
      local.signal <- getPhyloCmean(loc.tree.file, phenotype.table)
      stats.matrix[k.og, ] <- t(local.signal$stat)
      pval.matrix[k.og, ] <- t(local.signal$pvalue)
    }
    k.og <- k.og + 1
  }
  rownames(stats.matrix) <- og.names
  rownames(pval.matrix) <- og.names
  pval.matrix <- subset(pval.matrix, apply(X = pval.matrix, MARGIN = 1, FUN = sum)!=0)
  stats.matrix <- subset(stats.matrix, rownames(stats.matrix)%in%rownames(pval.matrix))
  ret.list <- list(stats = stats.matrix, pvals = pval.matrix)
  return(ret.list)
}
# 1.3 Initialize variables ===================================================#
# 1.3.1 Parser variables #
# Get script name
initial.options <- commandArgs(trailingOnly = FALSE)
script.name     <- basename(sub("--file=", "", initial.options[grep("--file=", initial.options)]))
# Process command line arguments
# Create a parser
option_list = list(
  make_option(c("-g", "--tree-dir"), type = "character", default = NULL,
              help ="Directory with the tree files for the OGs analyzed.", metavar = "character"),
  make_option(c("-p", "--phenotype_matrix"), type = "character", default = "NULL",
              help ="A file containing the phenotype values, each column is a different phenotype.", metavar = "character"),
  make_option(c("-s", "--min_num_strains"), type = "numeric", default = 30,
              help ="The minimum nuber of strains to include in the analysis, all trees with lower number of strains are dropped. [default= %default]", metavar = "numeric"),
  make_option(c("-o", "--output"), type = "character", default = "rank_list.txt",
              help ="Output filename [default= %default]", metavar = "character")
);
# Add command line arguments
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
opt_names <- names(opt)
if (is.null(opt$tree_dir)){
  print_help(opt_parser)
  err_str <- 'Argument missing "--tree_dir" must be provided.\n'
  stop(err_str, call.=FALSE)
} else if (is.null(opt$phenotype_matrix)){
  print_help(opt_parser)
  err_str <- 'Argument missing "--phenotype_matrix" must be provided.\n'
  stop(err_str, call.=FALSE)
}
# Parse the command line arguments
tree.dir <- opt$tree_dir
phenotype.matrix <- opt$phenotype_matrix
min.num.strains <- opt$min_num_strains
output <- opt$output
phenotype.table <- read.table(phenotype.matrix, sep = '\t', h = T, stringsAsFactors = F, row.names = 1)
#================== 2.0 Select the OGs by Cmean correlation ==================#
matrix.list <- getMatrices(tree.dir = tree.dir, phenotype.table = phenotype.table, tree.extension = '.out', min.num.strains = 30)
pvals <- matrix.list$pvals
stats <- matrix.list$stats
#=========================== 3.0 Write the output ============================#
write.table(x = pvals, file = 'pvals.out', sep = '\t', row.names = T, quote = F)
write.table(x = stats, file = 'stats.out', sep = '\t', row.names = T, quote = F)
