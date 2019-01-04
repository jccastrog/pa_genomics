/usr/bin/env Rscript

args <- commandArgs(TRUE)

ranks <- read.table(args[1], h = T, sep = '\t', stringsAsFactors = F)
rank.col <- ranks$rank
ranks <- data.matrix(ranks)
average_ranks <- function(rank.m){
	rank.col <- rank.m[,1]
	unique.ogs <- unique(sort(rank.m[,2]))
	num.ogs <- length(unique.ogs)
	num.folds <- dim(rank.m)[2] - 1
	ranks.mat <- matrix(NA, nrow = num.ogs, ncol = dim(rank.m)[2])
	ranks.mat[,1] <- unique.ogs
	for (i in 1:num.folds+1){
		for (j in 1:num.ogs){
			this.og <- rank.m[j,i]
			this.rank <- rank.col[j]
			ranks.mat[this.og,i] <- this.rank
		}
	}
	average.rank <- apply(ranks.mat, 1, mean)
	ranks.df <- data.frame(og = seq(1,num.ogs), rank = average.rank)
	ordered.ranks <- ranks.df[order(ranks.df$rank),]
	return(ordered.ranks)
}
	
rank.df <- average_ranks(ranks)

