#!/usr/bin/env Rscript

################################################################################
# Name:	    parse_SamplePanGenome.R
# Author:   Juan C. Castro [jcastro37@gatech.edu]
#           Georgia Institute of Technology
#           
# Version:  1.0
# Date:     05-Ago-2019
# License:  GNU General Public License v3.0.
# ==============================================================================
# 
################################################################################

#======================= 0.0 Install required packages ========================#
personal.lib.path = Sys.getenv("R_LIBS_USER")
if(!file.exists(personal.lib.path))
  dir.create(personal.lib.path)

packages <- c("gplots", "ggplot2", "optparse")
if(any(!(packages %in% installed.packages()))){
  cat("Please wait a moment! Installing required packages ...\n")
  install.packages(packages[!(packages %in% installed.packages())],
                   quiet = T, repos="http://cran.rstudio.com/",
                   lib = personal.lib.path)
  cat("Required packages installed!\n")
}
#======== 1.0 Load packages, define functions, and initialize variable========#
# 1.1 Load packages ==========================================================#
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(optparse))
# 1.2 Define functions =======================================================#
#' Take a pangenome matrix and convert it to numerical binary
#'
#' @param in.file Path to the file containing the matrix 
#' @return og.matrix A matrix with the information regarding the OGs in binary
#'                   terms where 1 is precense of the OG and 0 mean abscence 
#' @author Juan C. Castro \email{jcastro37@gatech.edu}
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
#' Take a pangenome matrix and convert it to numerical
#'
#' @param in.file Path to the file containing the matrix 
#' @return og.matrix A matrix with the information regarding the OGs in numerical
#'                   values
#' @author Juan C. Castro \email{jcastro37@gatech.edu}
parseOGNumerical <- function(in.file){
  text.mat <- read.table(in.file, header = T, stringsAsFactors = F)
  num.genomes <- ncol(text.mat)
  num.ogs <- nrow(text.mat)
  genome.list <- colnames(text.mat)
  numerical.mat <- matrix(0, nrow = num.ogs, ncol = num.genomes)
  for (i in 1:num.ogs){
    for (j in 1:num.genomes){
      if (text.mat[i,j] == "-"){
        numerical.mat[i,j] <- 0
      } else {
        ogs <- strsplit(text.mat[i,j], ",")[[1]]
        numerical.mat[i,j] <- length(ogs)
      }
    }
  }
  colnames(numerical.mat) <- genome.list
  return(numerical.mat)
}
#' Take a pangenome matrix and sample it to obtain the pan-core curve
#'
#' @param in.file Path to the file containing the matrix 
#' @param samp.perc.retain Number of genomes 
#' @param num.reps Number of repetitions in each sampling 
#' @return og.df A data.frame with the information regarding the OGs in numericals
#' @author Juan C. Castro \email{jcastro37@gatech.edu}
samplePanGenome <- function(in.file, samp.perc.retain, num.reps){
  binary.mat <-  parseOGBinary(in.file)
  num.genomes <- ncol(binary.mat)
  num.ogs <- nrow(binary.mat)
  og.df <- data.frame(genomes = seq(1, num.genomes), core_avg = NA, core_sd = NA, pan_avg = NA, pan_sd = NA)
  for (i in 1:num.genomes){
       core <- c()
       pan <- c()
       genome.pick <- seq(1, num.genomes)
       for (j in 1:num.reps){
         samp.genomes <- sample(genome.pick, i, replace = F)
         samp.mat <- binary.mat[, samp.genomes]
         if (i == 1){
           sum.samp <- samp.mat
         } else {
           sum.samp <- apply(X = samp.mat, MARGIN = 1, FUN = sum)
         }
         loc.core <- sum(sum.samp >= samp.perc.retain*i)
         loc.pan <- sum(sum.samp > 0)
         core <- c(core, loc.core)
         pan <- c(pan, loc.pan)
       }
       og.df$core_avg[i] <- mean(core)
       og.df$core_sd[i] <- sd(core)
       og.df$pan_avg[i] <- mean(pan)
       og.df$pan_sd[i] <- sd(pan)
  }
  return(og.df)
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
              help ="The OGs matrix that contains the gene names for each genome.", metavar = "character"),
  make_option(c("-b", "--convert_binary"), type = "logical", default = FALSE,
              help ="Convert the OGs matrix to a binary matrix.", metavar = "logical"),
  make_option(c("-n", "--convert_numerical"), type = "logical", default = FALSE,
              help ="Convert the OGs matrix to a numerical matrix where the number represents the number of copies of the gene.", metavar = "logical"),
  make_option(c("-s", "--sample_pangenome"), type = "logical", default = FALSE,
              help ="Sample the pangenome to construct a rarefraction matrix.", metavar = "logical"),
  make_option(c("-p", "--samp_perc_retain"), type = "numeric", default = NULL,
              help ="Fraction of the genomes necessary to consider a gene part of the pangenome.", metavar = "numeric"),
  make_option(c("-r", "--num_reps"), type = "numeric", default = NULL,
              help ="Number of repeated sampling of the pangeome in each iteration.", metavar = "numeric"),
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
convert.binary <- opt$convert_binary
convert.numerical <- opt$convert_numerical
sample.pangenome <- opt$sample_pangenome
samp.perc.retain <- opt$samp_perc_retain
num.reps <- opt$num_reps
output <- opt$output
#============== 2.0 Execute conversion procedures ==============#
# 2.1 Convert to binary ========================================#
if (convert.binary){
  binary.matrix <- parseOGBinary(in.file = pangenome.matrix)
  write.table(x = binary.matrix, file = output, quote = F, sep = '\t', row.names = F)
}
# 2.2 Convert to numerical =====================================#
if (convert.numerical){
  numerical.matrix <- parseOGNumerical(in.file = pangenome.matrix)
  write.table(x = numerical.matrix, file = output, quote = F, sep = '\t', row.names = F)
}
#============== 3.0 Sample the pangenome ==============#
if (sample.pangenome){
  raref.curve <- samplePanGenome(in.file = pangenome.matrix, samp.perc.retain = samp.perc.retain, num.reps = num.reps)
  newdf <- data.frame(genomes = c(seq(1,138),seq(1,138)), val = c(raref.curve$core_avg, raref.curve$pan_avg),
                      sd = c(raref.curve$core_sd, raref.curve$pan_sd), set = c(rep("Core", 138), rep("Pan", 138)))
  p <- ggplot(newdf, aes(colour=set)) + geom_point(aes( x = genomes, y = val)) +
    geom_errorbar(aes(x = genomes, ymin  = val - sd, ymax = val + sd)) +
    scale_colour_manual(values = c("#00FF00", "#0000FF")) +
    xlab("Number of genomes") + ylab("Number of OGs") +
    theme_classic() + theme(legend.position = "none")
  #write.table(x = raref.curve, file = output, quote = F, sep = '\t', row.names = F)
}

pdf("~/projects/pa_genomics/plots/pangenomePA.pdf", width = 10, height = 10)
p <- ggplot(newdf, aes(colour=set)) + geom_point(aes( x = genomes, y = val)) +
  geom_errorbar(aes(x = genomes, ymin  = val- sd, ymax = val + sd)) +
  scale_colour_manual(values = c("#D7B30B", "#008B85")) +
  xlab("Number of genomes") + ylab("Number of OGs") +
  theme_classic() + theme(legend.position = "none")
p
dev.off()
#======================================================#
