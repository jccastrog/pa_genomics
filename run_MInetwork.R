#!/usr/bin/env Rscript

################################################################################
# Name:     run_BGLR-RFE.R
# Author:   Juan C. Castro [jcastro37@gatech.edu]
#           Georgia Institute of Technology
#           
# Version:  1.0
# Date:     07-Aug-2016
# License:  GNU General Public License v3.0.
# ==============================================================================
# 
################################################################################

#======================= 0.0 Install required packages ========================#
personal.lib.path = Sys.getenv("R_LIBS_USER")
if(!file.exists(personal.lib.path))
	            dir.create(personal.lib.path)

packages <- c("entropy", "igraph")
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
suppressPackageStartupMessages(library(igraph))
# 1.2 Define functions =======================================================#
calculate_mi_graph <- function(initial.df){
	data.dim <- dim(initial.df)
	num.strains <- data.dim[1]
	num.ogs <- data.dim[2]
	og.names <- colnames(initial.df)
	graph.df <- data.frame(og1 = c(), og2 = c(), mi = c()) 
	for (i in 1:num.ogs){
		for (j in 1:num.ogs){
			if (i<j){
				x <- initial.df[,i]
				y <- initial.df[,j]
				joint.df  <- as.data.frame(cbind(x,y))
				lenX <- length(x)
				lenY <- length(y)
				sumX <- sum(x)
				sumY <- sum(y)
				pX <- c((lenX-sumX)/lenX, sumX/lenX)
				pY <- c((lenY-sumY)/lenY, sumY/lenY)
				joint.prob <- matrix(0,ncol = 2, nrow = 2)
				joint.prob[1,1] <- sum((joint.df$x==0)&(joint.df$y==0))
				joint.prob[1,2] <- sum((joint.df$x==1)&(joint.df$y==0))
				joint.prob[2,1] <- sum((joint.df$x==0)&(joint.df$y==1))
					joint.prob[2,2] <- sum((joint.df$x==1)&(joint.df$y==1))
				joint.prob <- joint.prob/lenX
				Hx <- entropy.plugin(pX)
				Hy <- entropy.plugin(pY)
				Hxy <- entropy.plugin(joint.prob)
				miXY <- Hx + Hy - Hxy
				temp.df <- data.frame(og1  = og.names[i], og2 = og.names[j], mi = miXY)
				graph.df <- rbind(graph.df, temp.df)
			}
		}
	}
	return(graph.df)
}


args <- commandArgs(TRUE)
input.file <- args[1]
output.file <- args[2]
initial.df <- read.table(input.file, header =  T, sep = '\t', stringsAsFactors = F, row.names = 1)
start.time <- Sys.time()
graph <- calculate_mi_graph(initial.df)
end.time <- Sys.time()
time.taken <- end.time - start.time
cat(paste(time.taken,'\n',sep=''))
write.table(x = graph, file = output.file, quote = F, sep ='\t', row.names = F)
