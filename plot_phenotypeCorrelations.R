#!/usr/bin/env Rscript

################################################################################
# Name:     plot_phenotypeCorrelations.R
# Author:   Juan C. Castro [jcastro37@gatech.edu]
#           Georgia Institute of Technology
#           
# Version:  1.0
# Date:     13-Sep-2016
# License:  GNU General Public License v3.0.
# ==============================================================================
# 
################################################################################

#======================= 0.0 Install required packages ========================#
personal.lib.path = Sys.getenv("R_LIBS_USER")
if(!file.exists(personal.lib.path))
            dir.create(personal.lib.path)

packages <- c("optparse","ggplot2")
if(any(!(packages %in% installed.packages()))){
        cat("Please wait a moment! Installing required packages ...\n")
        install.packages(packages[!(packages %in% installed.packages())],
                         quiet = T, repos="http://cran.rstudio.com/",
                         lib = personal.lib.path)
        cat("Required packages installed!\n")
}

#======= 1.0 Load packages, define functions, and initialize variable ========#
# 1.1 Load packages ==========================================================#
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(optparse))
# 1.2 Define functions =======================================================#
#' Plots two variables according to its properties (boxplot for factors and
#' continuous, colums for factors and factors, and scatter for both continuous)
#'
#' @param x A variable to be plotted (vector)
#' @param y A variable to be plotted (vector)
#' @param xname The name of the x variable (character)
#' @param yname The name of the y variable (character)
#' @param ... Fraphical parameters
#' @return p A ggplot object with the two variables
#' @author Juan C. Castro \email{jcastro37@gatech.edu}
plot2vars <- function(x, y, xname, yname){
	if (is.character(x) | is.factor(x)){
		if (is.character(y) | is.factor(y)){
			df <- data.frame(x = x, y = y)
			p <- ggplot(df, aes(x)) 
			p <- p + geom_bar(aes(fill = y))
			p <- p + xlab(xname) + ylab('Count')
			p <- p + theme_bw()
			barplot(df, xlab='Counts', main=paste('Distribution of ',xname,' and ',yname,sep=''))
		} else if (is.numeric(y)){
			df <- data.frame(x = x, y = y)
                        p <- ggplot(df, aes(x = x, y = y))
                        p <- p + geom_boxplot(aes(fill = 'grey55'))
                        p <- p + geom_jitter()
                        p <- p + xlab(xname) + ylab(yname)
                        p <- p + theme_bw()
			boxplot(y~x, col = 'grey55', xlab= xname, ylab = yname)
			aovar <- aov(y~x)
			pval <- summary(aovar)[[1]][["Pr(>F)"]][1]
			text(x = 2, y = max(y) ,labels = paste('P value=',pval,sep = ''))
		}
	} else if (is.numeric(x)){
		if (is.character(y) | is.factor(y)){
			df <- data.frame(x = x, y = y)
                        p <- ggplot(df, aes(x = y, y = x))
                        p <- p + geom_boxplot(aes(fill = 'grey55'))
                        p <- p + geom_jitter()
                        p <- p + xlab(xname) + ylab(yname)
                        p <- p + theme_bw() +  coord_flip()
			boxplot(x~y, col = 'grey55', xlab= xname, ylab = yname)
			aovar <- aov(x~y)
			pval <- summary(aovar)[[1]][["Pr(>F)"]][1]
			text(x = 2, y = max(x) ,labels = paste('Pvalue=',pval,sep = ''))
		} else if (is.numeric(y)){
			df <- data.frame(x = x, y = y)
                        p <- ggplot(df, aes(x = x, y = y))
			p <- p + geom_point(colour = 'grey55')
			p <- p + xlab(xname) + ylab(yname)
                        p <- p + theme_bw() +  coord_flip()
			plot(x,y, pch = 19, col = 'grey55', xlab= xname, ylab = yname)
			cval <- cor(x,y)                         
			ranx <- range(x)[2] - range(x)[1]
			rany <- range(y)[2] - range(y)[1]
			text(max(x)-(ranx/8), max(y)-(rany/8),labels = paste('R2=',cval,sep = ''))
		}
	}
	return(p)
}
# 1.3 Initialize variables ===================================================#
# 1.3.1 Parser variables
# Create a parser
option_list = list(
        make_option(c("-i", "--in_file"), type="character", default=NULL,
                help="The input table were colums are variables and rows are observations.", metavar="character"),
        make_option(c("-o", "--output"), type="character", default="output.pdf",
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
in.data <- read.table(opt$in_file, h = T , sep = ',', stringsAsFactors = F)
#2.2 Calculate the size of the table =========================================#
nvars <- ncol(in.data)
nobs <- nrow(in.data)
nams <- colnames(in.data)
#============================= 3.0 Plot the data =============================#
cat('Plotting data...\n')
pdf(opt$output)
for (i in 1:nvars){
	for (j in 1:nvars){
		if (j > i){
			x <- in.data[,i]
			y <- in.data[,j]
			xname <- nams[i]
			yname <- nams[j]
			plot2vars(x, y, xname, yname)
		}
	}
}
graphics.off()
#==============================================================================#
