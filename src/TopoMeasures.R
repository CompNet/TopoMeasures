######################################################################
# these functions are related to the processing of performance measures
# (in the context of community detection)
# source("CommunityDetection/performances/processing.R")
######################################################################
CommunityDetection_performances_processing_R <- TRUE

library("igraph")
library("flexclust")
library("mclust")

source("CommonTools/general/vectors.R")

source("CommonTools/properties/embeddedness.R")

source("CommunityDetection/partitions/files.R")
source("CommunityDetection/performances/general.R")




# processes the confusion matrix for two partitions, i.e. the matrix whose cell (i,j) contains
# the number of instances estimated to belong to community i when they actually belong to community j.
# TODO DONE
process.confusion.matrix <- function(part1, part2, remove.singles=FALSE)
{	#cat("part1:\n");print(part1)
	#cat("part2:\n");print(part2)
	
	# init
	num.row <- max(part1)
	num.col <- max(part2)
	elt.nbr <- num.row*num.col
	#cat("dimensions: ",num.row,"x",num.col,"=",elt.nbr,"\n",sep="")
	
	# TODO dirty workaround for when there're too many single communities,
	# resulting in a too large confusion matrix.
	# we simply merge as many single communities as needed to get
	# a matrix R can handle.
	if(elt.nbr > .Machine$integer.max || remove.singles)
	{	# select the input vector with highest community count
		part <- NULL
		if(num.row > num.col)
		{	part <- part1
			lim <- num.col
		}
		else
		{	part <- part2
			lim <- num.row
		}
		# process how many coms should be removed
		#excess <- elt.nbr - .Machine$integer.max
		# get the single coms
		t <- table(part)
		singles <- sort(which(t==1))
		if(length(singles)>0)
		{	#nbr.removed <- min(excess%/%lim+1,length(singles)-1)
			#last.single.index <- length(singles) - nbr.removed
			#last.single <- singles[last.single.index]
			#temp <- as.list(singles[(last.single.index+1):length(singles)])
			nbr.removed <- length(singles) - 1
			last.single <- singles[1]
			temp <- as.list(singles[2:length(singles)])
			sapply(temp,function(x)
					{	#cat(x,"/",max(singles),"\n")
						part[part==x] <<- last.single
					})
			# renumber communities
			old.coms <- sort(unique(part))
			indices <- as.list(1:length(old.coms))
			sapply(indices,function(index)
					{	#cat(index,"/",length(old.coms),": old com ",old.coms[index]," -> new com ",index,"\n",sep="")
						if(old.coms[index] != index)
							part[part==old.coms[index]] <<- index
					})
		}
		if(num.row > num.col)
			part1 <- part
		else
			part2 <- part
		
		num.row <- max(part1)
		num.col <- max(part2)
		elt.nbr <- num.row*num.col
		#cat("new dimensions: ",num.row,"x",num.col,"=",elt.nbr,"\n",sep="")
	}
	
	# process matrix
#	conf.matrix <- matrix(0,nrow=num.row,ncol=num.col)
#	for(i in 1:length(part1))
#	{	row <- part1[i]
#		col <- part2[i]
#		conf.matrix[row,col] <- conf.matrix[row,col] + 1
#	}
	
	
	# more efficient:
	conf.matrix <- as.matrix(table(part1,part2)) # note: some rows/cols can be missing if they don't appear in one partition
	
	return(conf.matrix);
}



# processes the pairwise matrix for two partitions, i.e. the matrix whose cells contain:
#	- (0,0): number of pairs of instances classified in the same community by both methods
#	- (1,1): number of pairs of instances classified in different communities by both methods
#	- (0,1),(1.0): number of pairs of instances classified in different communities by one method, and in the same by the other
# TODO DONE
process.pairwise.matrix <- function(part1, part2)
{	# init
	result <- matrix(0,nrow=2,ncol=2)
	
	# process matrix
	for(i in 1:(length(part1)-1))
	{	for(j in (i+1):length(part1))
		{	#cat("i,j:",i,",",j," 1:",part1[i],",",part1[j]," 2:",part2[i],",",part2[j],"\n")			
			if(part1[i]==part1[j])
			{	if(part2[i]==part2[j])
				{	# both the same
					#cat("both same\n")				
					result[1,1] <- result[1,1] + 1
				}
				else
				{	# the same for est., different for ref.
					#cat("same diff\n")				
					result[1,2] <- result[1,2] + 1
				}
			}
			else
			{	if(part2[i]==part2[j])
				{	# different for est., the same for ref.
					#cat("diff same\n")				
					result[2,1] <- result[2,1] + 1
				}
				else
				{	# different for both
					#cat("both diff\n")				
					result[2,2] <- result[2,2] + 1
				}
			}
		}
	}
	
	return(result);
}



# processes the percent of correctly classified nodes for the specified partitions
# it identifies similar communities by selecting the max confusion matrix cell (i,j):
# then estimated community i is supposed to correspond to actual community class j.
# this may not be the best way...
# TODO DONE
process.percent.correct<-function(reference, estimation, remove.singles=FALSE)
{	# process the confusion matrix
	conf.matrix <- process.confusion.matrix(reference,estimation,remove.singles);
	#print(conf.matrix)
	
	norm <- sum(conf.matrix)
	inter.nbr <- min(dim(conf.matrix))
	#cat("norm:",norm,"\n")
	#cat("inter.nbr:",inter.nbr,"\n")
	#cat("dim:",dim(conf.matrix),"\n")
	
	# match the communities
	ref.used <- c();
	est.used <- c();
	total <- 0;
	val <- 1
	i <- 1
	while(i<=inter.nbr && val>0)
	{	#cat("percent.correct: ",i,"/",inter.nbr,"\n",sep="")
		mx <- get.max.mat(conf.matrix)
		val <- conf.matrix[mx[1],mx[2]]
		#cat("i:",i," val:",val," mx:",mx,"\n")
		if(val>0)
		{	if(is.na(get.last.position(mx[1],ref.used)) && is.na(get.last.position(mx[2],est.used)))
			{	ref.used <- c(ref.used,mx[1])
				est.used <- c(est.used,mx[2])
				total <- total + conf.matrix[mx[1],mx[2]]
			}
			conf.matrix[mx[1],mx[2]] <- 0
			i <- i + 1
		}
	}
	
	# process final result
	result <- total/norm;
	return(result);
}



# like process.percent.correct, except it includes the
# variant by newman, in which if two reference communities
# correspond to the same estimated community, then all nodes
# are considered as misclassified (option no.merge).
# Newman uses the inverted version of the measure, i.e. he
# considers the purity of real communities (instead of the estimated ones)
process.purity <- function(reference, estimation, inverted=FALSE, remove.singles=FALSE, no.merge=FALSE)
{	# possibly invert the vectors
	if(inverted)
	{	temp <- reference
		reference <- estimation
		estimation <- temp
	}	
	
	# process the confusion matrix
	conf.matrix <- process.confusion.matrix(estimation,reference,remove.singles)
	#print(dim(conf.matrix))
	#print(conf.matrix)
	
	# init
	norm <- sum(conf.matrix)
	corresp <- rep(0,nrow(conf.matrix))
	#print(conf.matrix)
	
	# for each estimated community, identify the corresponding reference community
	for(e in 1:nrow(conf.matrix))
		corresp[e] <- which.max(conf.matrix[e,])
	#print(corresp)
	
	# identify and count the correctly classified nodes
	correct <- 0
	for(r in 1:ncol(conf.matrix))
	{	indices <- which(corresp==r)
		#print(indices)
		# possibly ignore merged communities
		if(length(indices)==1 || !no.merge)
			correct <- correct + sum(conf.matrix[indices,r])
		#cat("correct:",correct,"\n")
	}
	
	result <- correct / norm
	return(result)
}



# variant of the purity taking topological properties
# into account. the parameter topo.measure associates
# some topological property to each node. this information
# is used to determine how much the classifier must be
# penalized in case of error.
process.topological.purity <- function(reference, estimation, inverted=FALSE, topo.measure, remove.singles=FALSE, details=FALSE)
{	# possibly invert the vectors
	if(inverted)
	{	temp <- reference
		reference <- estimation
		estimation <- temp
	}
	
	# process the confusion matrix
	conf.matrix <- process.confusion.matrix(estimation,reference,remove.singles)
	#print(dim(conf.matrix))
	
	# init
	norm <- sum(topo.measure)
	topo.measure <- topo.measure / norm
	#cat("norm:",norm,"\n",ep="")
	corresp <- rep(0,nrow(conf.matrix))
	
	# for each reference community, identify the corresponding estimated community
	for(e in 1:nrow(conf.matrix))
		corresp[e] <- which.max(conf.matrix[e,])
	#print(corresp)
	
	# identify and count the correctly classified nodes
	correct <- 0
	correctness <- rep(NA,length(reference))
	for(i in 1:length(reference))
	{	#cat("estimation:",estimation[i],"\n")
		#cat("reference:",reference[i],"\n")
		corres <- corresp[estimation[i]]
		#cat("corres:",corres,"\n")
		if(!is.na(corres) && reference[i]==corres)
		{	correct <- correct + topo.measure[i]
			#cat("correct:",topo.measure[i],"\n",sep="")
			correctness[i] <- 1
		}
		else
		{	#cat("INcorrect:",topo.measure[i],"\n",sep="")
			correctness[i] <- 0
		}
	}
	
	if(details)
		result <- correctness
	else
		result <- correct
	#cat("correct.old:",correct.old/length(estimation),"\n")
	#cat("purity",process.purity(reference, estimation, inverted),"\n")
	#cat("percent.correct",process.percent.correct(reference,estimation),"\n")
	#cat("result:",result,"\n")
	#cat("xxxxxxxxxxxxxxxxxxxxxx\n")
	#print(correctness)
	return(result)
	
# tests	
#	conf.matrix <- matrix(c(12,2,3,5,16,8,9,3,20,5,6,18),ncol=4)
#	reference <- c(1,1,1,1,1,1,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4)
#	estimation <- c(1,1,1,1,2,3,1,2,2,2,2,2,3,1,2,3,3,3,3,3,3,1,2,3,3,3,3,3,3,3)
}




# processes the normalized mutual information measure for the specified partitions
# TODO DONE
process.NMI <- function(part1, part2, remove.singles=FALSE)
{	# process the confusion matrix
	conf.matrix <- process.confusion.matrix(part1,part2,remove.singles)
	num.row <- nrow(conf.matrix)
	num.col <- ncol(conf.matrix)
	
	# process all sums once and for all
	n <- sum(conf.matrix)
	s.row <- rowSums(conf.matrix)
	s.col <- colSums(conf.matrix)
	
	# process normalized mutual information measure
	sum1 <- 0;
	for(i in 1:num.row)
	{	for(j in 1:num.col)
		{	temp <- conf.matrix[i,j] * log2((conf.matrix[i,j]*n) / (s.row[i]*s.col[j]))
			if(!is.na(temp))
				sum1 <- sum1 + temp;
		}
	}
	
	sum2 <- 0;
	for(i in 1:num.row)
	{	temp <- s.row[i]*log2(s.row[i]/n)
		if(!is.na(temp))
			sum2 <- sum2 + temp;
	}
	
	sum3 <- 0;
	for(j in 1:num.col)
	{	temp <- s.col[j]*log2(s.col[j]/n)
		if(!is.na(temp))
			sum3 <- sum3 + temp;
	}
	
	# combine all the parts...
	info.metric <- (-2)*sum1 / (sum2+sum3)
	return(info.metric)
}



process.topological.NMI <- function(part1, part2, topo.measure, remove.singles=FALSE)
{	norm <- sum(topo.measure)
	topo.measure <- topo.measure/norm
	
	# process the (topological-weighted) confusion matrix
	num.row <- max(part1)
	num.col <- max(part2)
	conf.matrix <- matrix(0,nrow=num.row,ncol=num.col)
	for(i in 1:length(part1))
	{	row <- part1[i]
		col <- part2[i]
		value <- topo.measure[i]
		conf.matrix[row,col] <- conf.matrix[row,col] + value
	}
	#print(conf.matrix)	
	
	# process all sums once and for all
	s.row <- rowSums(conf.matrix)
	s.col <- colSums(conf.matrix)
	
	# process normalized mutual information measure
	sum1 <- 0;
	for(i in 1:num.row)
	{	for(j in 1:num.col)
		{	temp <- conf.matrix[i,j] * log2((conf.matrix[i,j]) / (s.row[i]*s.col[j]))
			if(!is.na(temp))
				sum1 <- sum1 + temp;
		}
	}
	
	sum2 <- 0;
	for(i in 1:num.row)
	{	temp <- s.row[i]*log2(s.row[i])
		if(!is.na(temp))
			sum2 <- sum2 + temp;
	}
	
	sum3 <- 0;
	for(j in 1:num.col)
	{	temp <- s.col[j]*log2(s.col[j])
		if(!is.na(temp))
			sum3 <- sum3 + temp;
	}
	
	# combine all the parts...
	info.metric <- (-2)*sum1 / (sum2+sum3)
	return(info.metric)
}



# processes the rand index for the specified partitions
# the rand index is nothing but the percent correct 
# processed on the pairwise matrix
# TODO DONE
process.rand.index<-function(part1, part2)
{	#mat <- process.pairwise.matrix(part1,part2)
	#result <- (mat[1,1]+mat[2,2]) / sum(mat)
	#cat(paste("RI: original=",result,sep=""))
	
	# more efficient
	#conf.matrix <- process.confusion.matrix(part1,part2)
	#result <- randIndex(conf.matrix,correct=FALSE)
	result <- randIndex(as.integer(part1),as.integer(part2),correct=FALSE)
	#cat(paste("new=",result,"\n",sep=""))
	
	return(result)
}



process.topological.rand.index<-function(part1, part2, topo.measure)
{	topo.measure <- topo.measure / sum(topo.measure)
	mat <- matrix(0,nrow=2,ncol=2)
	
	# process matrix
	for(i in 1:(length(part1)-1))
	{	for(j in (i+1):length(part1))
		{	#cat("i,j:",i,",",j," 1:",part1[i],",",part1[j]," 2:",part2[i],",",part2[j],"\n")			
			value <- topo.measure[i] * topo.measure[j]
			if(part1[i]==part1[j])
			{	if(part2[i]==part2[j])
				{	# both the same
					#cat("both same\n")				
					mat[1,1] <- mat[1,1] + value
				}
				else
				{	# the same for est., different for ref.
					#cat("same diff\n")				
					mat[1,2] <- mat[1,2] + value
				}
			}
			else
			{	if(part2[i]==part2[j])
				{	# different for est., the same for ref.
					#cat("diff same\n")				
					mat[2,1] <- mat[2,1] + value
				}
				else
				{	# different for both
					#cat("both diff\n")				
					mat[2,2] <- mat[2,2] + value
				}
			}
		}
	}
	#print(mat)
	
	result <- (mat[1,1]+mat[2,2]) / sum(mat)
	#cat(paste("RI: original=",result,sep=""))
	
	return(result)
}




# processes the adjusted rand index for the specified partitions
# the rand index is nothing but the kappa statistic 
# processed on the pairwise matrix
# TODO DONE
process.adjusted.rand.index<-function(part1, part2)
{	#result <- adjustedRandIndex(part1,part2)
	#cat(paste("ARI: original=",result,sep=""))
	
	# alternative
	#conf.matrix <- process.confusion.matrix(part1,part2);
	#result <- randIndex(conf.matrix,correct=TRUE)
	result <- randIndex(as.integer(part1),as.integer(part2),correct=TRUE)
	#cat(paste("new=",result,"\n",sep=""))
	
	return(result)
}




process.topological.adjusted.rand.index<-function(part1, part2, topo.measure)
{	# whole graph
	ws <- generalized.coefficient(1:length(part1),topo.measure)
	
	# parts 1
	s1 <- max(part1)
	w1 <- 0
	for(i in 1:s1)
	{	indices1 <- which(part1==i)
		w1 <- w1 + generalized.coefficient(indices1,topo.measure)
	}
	
	# parts 2
	s2 <- max(part2)
	w2 <- 0
	for(j in 1:s2)
	{	indices2 <- which(part2==j)
		w2 <- w2 + generalized.coefficient(indices2,topo.measure)
	}
	
	# intersections
	w12 <- 0
	for(i in 1:s1)
	{	indices1 <- which(part1==i)
		for(j in 1:s2)
		{	indices2 <- which(part2==j)
			indices12 <- intersect(indices1,indices2)
			w12 <- w12 + generalized.coefficient(indices12,topo.measure)
		}
	}
	
	#cat("ws=",ws," w1=",w1," w2=",w2," w12=",w12,"\n",sep="")
	
	# process measure (cf. the formula in the paper)
	result <- (w12 - w1*w2/ws) / (1/2*(w1+w2) - w1*w2/ws)
	return(result)
}




# process the generalized binomial coefficient used in the adjusted
# rand index calculation.
# the indices are positions in the whole group of nodes
generalized.coefficient <- function(indices, topo.measure)
{	# init
	size <- length(indices)
	result <- 0
	#topo.measure <- rep(1,length(topo.measure)) #check comparison with regular ARI
	
	if(size>1)
	{	# consider each pair and their weight
		for(i in 1:(size-1))
		{	for(j in (i+1):size)
			{	idx1 <- indices[i]
				idx2 <- indices[j]
				result <- result + topo.measure[idx1]*topo.measure[idx2]
				#cat(topo.measure[idx1]*topo.measure[idx2],":",sep="")
			}
		}
	}
	#cat("\n")
	
	#print(indices)
	#print(topo.measure)
	#print(result)
	#print("-------------------------------------")
	
	return(result)
}




# determines the partition corresponding to the best cut.
# two different criteria can be used:
#	- maximal modularity first:	the partition with the highest modularity is returned
# 	- smaller partition first: 	begining with the partition with the fewest communities,
#								the function goes up in the hierarchy and stops when the modularity
#								does not improve significantly. It means the partition is not necessarily
#								the best in terms of modularity, but it still has a significantly similar
#								modularity and less communities (which is usually a good thing)
cut.dendrogram <- function(g, result, algoname, criterion=c("modularity","size"), verbose=0, logfile)
{	# possibly process the modularities (when the algo didn't)
	size <- dim(result$merges)[1];
	if(is.null(result$modularity))
	{	result$modularity <- rep(0,size+1);
		# special case: flat dendrogram 
		if(size==0)
		{	membership <- result$membership;
			if(verbose>0 && !is.null(logfile))
				cat("..WARNING: empty dendrogram!\n",file=logfile,append=TRUE,sep="");
		}
		# normal case
		else
		{	for(i in 0:size)
			{	if(algoname=="EIGENVECTOR")
					membership <- community.le.to.membership(merges=result$merges,steps=i,membership=result$membership)$membership
				else
					membership <- community.to.membership(g,merges=result$merges,steps=i)$membership
				result$modularity[i+1] <- modularity(g,membership);
			}
		}
	}
	
	# get the best step
	if(criterion=="modularity")
	{	# max modularity
		best.idx <- which.max(result$modularity)
	}
	else if(criterion=="size")
	{	# min size
		best.idx <- get.smaller.similar.value(values=result$modularity)
	}
	steps <- best.idx - 1
	
	# get the corresponding partition
	if(steps>0)
	{	if(algoname=="EIGENVECTOR")
			membership <- community.le.to.membership(merges=result$merges,steps=steps,membership=result$membership)$membership
		else
			membership <- community.to.membership(g,merges=result$merges,steps=steps,membership=TRUE,csize=FALSE)$membership
	}
	
	return(membership);
}



# return the partition corresponding to the specified level in the hierarchy
# represented by the specified dendrogram.
# note: the first level is noted 0.
get.level.partition <- function(g, algoname, dendrogram, level=0)
{	
	size <- dim(dendrogram$merges)[1]
	
	# level depeer than the dendrogram height (ERROR)
	if(level>size)
	{	result <- NULL
	}
	
	# appropriate level
	else
	{	# special case: flat dendrogram 
		if(size==0)
		{	result <- dendrogram$membership
		}
		
		# special case: level 0 (no merge)
		else if(level==0)
		{	if(algoname=="EIGENVECTOR")
				result <- community.le.to.membership(merges=dendrogram$merges,steps=level,membership=dendrogram$membership)$membership
			else
				result <- community.to.membership(g,merges=dendrogram$merges,steps=level)$membership
		}
		
		# normal case: some other level requiring to merge
		else
		{	if(algoname=="EIGENVECTOR")
			{	result <- community.le.to.membership(merges=dendrogram$merges,steps=level,membership=dendrogram$membership)$membership
			}
			else
			{	# there seems to be a bug with walktrap, which sometimes generates dendrograms with (0 0) lines,
				# meaning impossible merge (the dendrogram being here a matrix with two columns corresponding to the 
				# various merging steps, each couple of values representing the communities to be merged)
				if(!(dendrogram$merges[level,1]<1 && dendrogram$merges[level,2]<1))
					result <- community.to.membership(g,merges=dendrogram$merges,steps=level)$membership
				else
				{	cat("WARNING: the algorithm (certainly Walktrap) generated a dendrogram with a (0 0) line, which it should not!\n")
					result <- dendrogram$membership
				}
			}
		}
	}
	
	return(result)
}



# return all the partitions from a hierarchy, under the form of a list of membership vectors
# TODO DONE
get.all.partitions <- function(g, algoname, dendrogram)
{	result <- list()
	
#print(g)
#print(algoname)
#print(dendrogram)
	
	size <- dim(dendrogram$merges)[1]
#cat("1\n")	
	# special case: flat dendrogram 
	if(size==0)
	{	
#cat("2\n")
		membership <- dendrogram$membership
#cat("3\n")
		result[[1]] <- membership
#cat("4\n")
	}
	# normal case
	else
	{	
#cat("5\n")	
		for(i in 0:size)
		{	
#cat("6:",i,"\n")	
			if(algoname=="EIGENVECTOR")
			{	
#cat("7:",i,"\n")	
				membership <- community.le.to.membership(merges=dendrogram$merges,steps=i,membership=dendrogram$membership)$membership
#cat("8:",i,"\n")	
			}
			
			# special case: level 0 (no merge)
			else if(i==0)
			{	if(algoname=="EIGENVECTOR")
					membership <- community.le.to.membership(merges=dendrogram$merges,steps=i,membership=dendrogram$membership)$membership
				else
					membership <- community.to.membership(g,merges=dendrogram$merges,steps=i)$membership
			}
			
			# normal case: some other level requiring to merge
			else
			{	
#cat("9:",i,"\n")
				if(!(dendrogram$merges[i,1]<1 && dendrogram$merges[i,2]<1))
					membership <- community.to.membership(g,merges=dendrogram$merges,steps=i)$membership
				else
				{	cat("WARNING: the algorithm (certainly Walktrap) generated a dendrogram with a (0 0) line, which it should not!\n")
					membership <- dendrogram$membership
				}
#cat("10:",i,"\n")
			}
#cat("11:",i,"\n")
#print(membership)
			result[[as.character(i)]] <- membership
#cat("12:",i,"\n")
#print(dim(dendrogram$merges))
#print(dendrogram$merges)
		}
#cat("13\n")
	}
#cat("14\n")
	
	return(result)
}



# processes the performance for the specified metric
# TODO DONE
process.performance.metric <- function(g=NA, partition, reference, perfname, remove.singles=FALSE)
{	result <- NA
	
	if(perfname=="MODULARITY")
	{	#if(!is.na(g))
		result <- modularity(g,partition)
		#else
		#	result <- NA
	}
	else if(perfname=="NORM_MUT_INF" && !is.na(reference))
		result <- process.NMI(partition,reference,remove.singles)
	
	else if(perfname=="PRCNT_CRCT" && !is.na(reference))
		result <- process.percent.correct(reference, estimation=partition, remove.singles)
	else if(perfname=="PURITY_REF" && !is.na(reference))
		result <- process.purity(reference, estimation=partition, inverted=TRUE, remove.singles)
	else if(perfname=="PURITY_EST" && !is.na(reference))
		result <- process.purity(reference, estimation=partition, inverted=FALSE, remove.singles)
	else if(perfname=="PURITY_BOTH")
	{	result1 <- process.purity(reference, estimation=partition, inverted=FALSE, remove.singles)
		result2 <- process.purity(reference, estimation=partition, inverted=TRUE, remove.singles)
		alpha <- 0.5
		result <- 1 / (alpha/result1 + (1-alpha)/result2)
	}
	
	else if(perfname=="RAND_IDX" && !is.na(reference))
		result <- process.rand.index(partition,reference)
	else if(perfname=="ADJ_RAND_IDX" && !is.na(reference))
		result <- process.adjusted.rand.index(partition,reference)
	
	else if(perfname=="TOPO_PURITY_REF_DEGREE")
		result <- process.topological.purity(reference, estimation=partition, 
				topo.measure=degree(g)/max(degree(g)), inverted=TRUE, remove.singles)
	else if(perfname=="TOPO_PURITY_REF_EMBEDDEDNESS")
		result <- process.topological.purity(reference, estimation=partition, 
				topo.measure=embeddedness(g,reference), inverted=TRUE, remove.singles)
	else if(perfname=="TOPO_PURITY_REF_BOTH")
		result <- process.topological.purity(reference, estimation=partition, 
				topo.measure=degree(g)/max(degree(g))*embeddedness(g,reference), inverted=TRUE, remove.singles)
	else if(perfname=="TOPO_PURITY_EST_DEGREE")
		result <- process.topological.purity(reference, estimation=partition, 
				topo.measure=degree(g)/max(degree(g)), inverted=FALSE, remove.singles)
	else if(perfname=="TOPO_PURITY_EST_EMBEDDEDNESS")
		result <- process.topological.purity(reference, estimation=partition, 
				topo.measure=embeddedness(g,reference), inverted=FALSE, remove.singles)
	else if(perfname=="TOPO_PURITY_EST_BOTH")
		result <- process.topological.purity(reference, estimation=partition, 
				topo.measure=degree(g)/max(degree(g))*embeddedness(g,reference), inverted=FALSE, remove.singles)
	else if(perfname=="TOPO_PURITY_BOTH_DEGREE")
	{	result1 <- process.topological.purity(reference, estimation=partition, 
				topo.measure=degree(g)/max(degree(g)), inverted=FALSE, remove.singles)
		result2 <- process.topological.purity(reference, estimation=partition, 
				topo.measure=degree(g)/max(degree(g)), inverted=TRUE, remove.singles)
		alpha <- 0.5
		result <- 1 / (alpha/result1 + (1-alpha)/result2)
	}
	else if(perfname=="TOPO_PURITY_BOTH_EMBEDDEDNESS")
	{	result1 <- process.topological.purity(reference, estimation=partition,
				topo.measure=embeddedness(g,reference), inverted=FALSE, remove.singles)
		result2 <- process.topological.purity(reference, estimation=partition,
				topo.measure=embeddedness(g,reference), inverted=TRUE, remove.singles)
		alpha <- 0.5
		result <- 1 / (alpha/result1 + (1-alpha)/result2)
	}
	else if(perfname=="TOPO_PURITY_BOTH_BOTH")
	{	result1 <- process.topological.purity(reference, estimation=partition, 
				topo.measure=degree(g)/max(degree(g))*embeddedness(g,reference), inverted=FALSE, remove.singles)
		result2 <- process.topological.purity(reference, estimation=partition, 
				topo.measure=degree(g)/max(degree(g))*embeddedness(g,reference), inverted=TRUE, remove.singles)
		alpha <- 0.5
		result <- 1 / (alpha/result1 + (1-alpha)/result2)
	}
	else if(perfname=="TOPO_PURITY_BOTH_CLOSENESS")
	{	result1 <- process.topological.purity(reference, estimation=partition,
				topo.measure=closeness(g), inverted=FALSE, remove.singles)
		result2 <- process.topological.purity(reference, estimation=partition,
				topo.measure=closeness(g), inverted=TRUE, remove.singles)
		alpha <- 0.5
		result <- 1 / (alpha/result1 + (1-alpha)/result2)
	}
	
	
	else if(perfname=="TOPO_ADJ_RAND_IDX_DEGREE")
		result <- process.topological.adjusted.rand.index(part1=reference, part2=partition, 
				topo.measure=degree(g)/max(degree(g)))
	else if(perfname=="TOPO_ADJ_RAND_IDX_EMBEDDEDNESS")
		result <- process.topological.adjusted.rand.index(part1=reference, part2=partition, 
				topo.measure=embeddedness(g,reference))
	else if(perfname=="TOPO_ADJ_RAND_IDX_BOTH")
		result <- process.topological.adjusted.rand.index(part1=reference, part2=partition, 
				topo.measure=degree(g)/max(degree(g))*embeddedness(g,reference))
	else if(perfname=="TOPO_RAND_IDX_DEGREE")
		result <- process.topological.rand.index(part1=reference, part2=partition, 
				topo.measure=degree(g)/max(degree(g)))
	else if(perfname=="TOPO_RAND_IDX_EMBEDDEDNESS")
		result <- process.topological.rand.index(part1=reference, part2=partition, 
				topo.measure=embeddedness(g,reference))
	else if(perfname=="TOPO_RAND_IDX_BOTH")
		result <- process.topological.rand.index(part1=reference, part2=partition, 
				topo.measure=degree(g)/max(degree(g))*embeddedness(g,reference))
	
	else if(perfname=="TOPO_NORM_MUT_INF_DEGREE")
		result <- process.topological.NMI(part1=reference, part2=partition, 
				topo.measure=degree(g)/max(degree(g)), remove.singles)
	else if(perfname=="TOPO_NORM_MUT_INF_EMBEDDEDNESS")
		result <- process.topological.NMI(part1=reference, part2=partition, 
				topo.measure=embeddedness(g,reference), remove.singles)
	else if(perfname=="TOPO_NORM_MUT_INF_BOTH")
		result <- process.topological.NMI(part1=reference, part2=partition, 
				topo.measure=degree(g)/max(degree(g))*embeddedness(g,reference), remove.singles)
	
	return(result);
}



# processes all performance measures from a network, a partition and a reference
process.performance <- function(g, filename1, filename2, performances.names=c(), remove.singles=FALSE)
{	# all performances if performances.names is empty
	if(length(performances.names)==0)
		performances.names <- get.performances.names()
	
	# read the partition files
	#cat("estimation:",filename1,"\n")	
	#cat("reference:",filename2,"\n")	
	partition <- load.partition(filename1)
	if(is.na(filename2))
		reference <- NA
	else
		reference <- load.partition(filename2)
	
	# process the performances
	result <- rep(NA,length(performances.names))
	for(i in 1:length(performances.names))
	{	#cat("i:",i," perf:",performances.names[i],"\n",sep="")
		result[i] <- process.performance.metric(g,partition,reference,performances.names[i],remove.singles)
	}
	
	return(result)
}



# process performances for only some of the actual communities
process.restricted.performance <- function(part.filename, ref.filename, performances.names=c(), bins=c())
{	# init performance names
	if(length(performances.names)==0)
		performances.names <- get.performances.names();
	
	# init or normalize bins
	if(length(bins)==0)
		bins <- c(0.4,0.3,0.3)
	else
		bins <- bins/sum(bins)
	
	# read the partition files
	partition <- load.partition(part.filename);
	reference <- load.partition(ref.filename);
	#partition <- c(1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,4,4,4,5,5,6,6)
	#reference <- c(1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,3,3,4,4,5,5)
	
	# process community sizes
	n <- length(reference)
	com.count <- max(reference)
	com.freq <- rep(0,com.count)
	for(i in 1:com.count)
		com.freq[i] <- sum(reference==i)
	#cat(com.count,"\n")
	#cat(com.freq,"\n")
	
	
	# init sub-partitions
	indices <- list()
	temp.n <- 0
	temp.b <- 0
	freq <- 1
	for(i in 1:length(bins))
	{	temp.idx <- c()
		temp.b <- temp.b + bins[i]
		while(temp.n<(temp.b*n) && freq>0)
		{	# get the biggest remaining com
			max.com <- which.max(com.freq)
			# grab its size
			freq <- com.freq[max.com]
			if(freq>0)
			{	# add to the current total
				temp.n <- temp.n + freq
				# reset the com size
				com.freq[max.com] <- 0
				# get the nodes in this com
				nodes <- which(reference==max.com)
				# add them to the current subpartition
				temp.idx <- c(temp.idx,nodes)
			}
		}
		indices[[i]] <- temp.idx
	}
	#print(indices)
	
	# process restricted performances
	result <- list(bin=round(100*bins),perf=matrix(NA,nrow=length(bins),ncol=length(performances.names)))
	for(i in 1:length(bins))
	{	if(length(indices[[i]]>0))
		{	subpartition <- partition[indices[[i]]]
			#print(subpartition)			
			subreference <- reference[indices[[i]]]
			#print(subreference)			
			for(j in 1:length(performances.names))
				result$perf[i,j] <- process.performance.metric(NA,subpartition,subreference,performances.names[j]);
		}
	}
	
	return(result);
}
