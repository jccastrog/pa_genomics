#!/usr/bin/env Rscript

################################################################################
# Name:     parse_prepareData.R
# Author:   Juan C. Castro [jcastro37@gatech.edu]
#           Georgia Institute of Technology
#           
# Version:  1.0
# Date:     04-Dic-2017
# License:  GNU General Public License v3.0.
# ==============================================================================
# 
################################################################################

#======================= 0.0 Install required packages ========================#
personal.lib.path = Sys.getenv("R_LIBS_USER")
if(!file.exists(personal.lib.path))
	            dir.create(personal.lib.path)

packages <- c("optparse")
if(any(!(packages %in% installed.packages()))){
	cat("Please wait a moment! Installing required packages ...\n")
        install.packages(packages[!(packages %in% installed.packages())],
			 quiet = T, repos="http://cran.rstudio.com/",
			 lib = personal.lib.path)
	cat("Required packages installed!\n")
}

#======== 1.0 Load packages, define functions, and initialize variable========#
# 1.1 Load packages ==========================================================#
library(optparse)
# 1.2 Define functions =======================================================#
#' Replaces the character values for 1s and NA valuess foe 0s
#'
#' @param data A matrix with character values to be replaced
#' @return transform A matrix with 1 and 0 values
#' @author Juan C. Castro \email{jcastro37@gatech.edu}
replaceString2Num <- function(data){	
	mat.data <- as.matrix(data)
	mat.data[mat.data == ""] <- NA
	mat.data[mat.data == "-"] <- NA
	mat.data[is.na(mat.data)] <- 0
	mat.data[mat.data != 0] <- 1
	transform <- apply(mat.data, 2, as.numeric)
	return(transform)
}
# 1.3 Initialize variables ===================================================#
# 1.3.1 Parser variables #
# Create a parser
option_list = list(
		   make_option(c("-i", "--in_file"), type="character", default=NULL,
			       help="The input table were colums are variables and rows are observations.", metavar="character"),
		   make_option(c("-o", "--output"), type="character", default="output.tsv",
			       help="Output file name [default= %default]", metavar="character")
		   );
# Parse the command line arguments
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
opt_names <- names(opt)
if (is.null(opt$in_file)){
	print_help(opt_parser)
	err_str <- 'Argument missing "--in_file" must be provided.\n'
	stop(err_str, call.=FALSE)
}

#======================== 2.0 Load and parse the data ========================#
#======================== 2.0 Load and parse the data ========================#
# 2.1 Load the data ==========================================================#
in.data <- read.table(opt$in_file, h = T , sep = '\t', stringsAsFactors = F)
cnams <- colnames(in.data)
#=========================== 3.0 Transform the data ==========================#
# 3.1 Apply the transformation function ======================================#
out.data <- replaceString2Num(in.data)
#colnames(out.data) <- cnams
# 3.2 Write the new data =====================================================#
write.table(x = out.data, file = opt$output, sep = '\t', quote = F, row.names = F)
#=============================================================================#
