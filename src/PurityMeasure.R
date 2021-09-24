############################################################################
# Functions used to process the purity measure and its variants.
#
# Vincent Labatut 2012-16
############################################################################
# get the common functions
source("src/CommonFunctions.R")




############################################################################
# Processes the proportion of similarly classified elements for the specified 
# partitions.
# 
# In order to identify similar parts in two distinct partitions, the function
# selects the cell (i,j) of the confusion matrix which has the maximal value:
# part #i in partition #1 is supposed to correspond to part #j in partition #2.
# Then, all remaining values in row #i and column #j are discarded (since the
# association has been set, and we look for a 1-to-1 mapping). We go one similarly
# with the rest of the matrix, until there is no more value to process. If the
# numbers of parts are different in the two considered partitions, the smallest
# parts of the larger partition will be left over (i.e. no mapping for them).
#
# This measure is found quite commonly in the cluster analysis literature, and
# in the early community detection works. 
#
# partition1: the first partition to consider, represented as an integer vector.
#			  Each value in the vector represents the id of a part. The parts 
#			  must be counted starting from one (not zero).
# partition2: the second partition to consider. Same representation as for
#			  the first one. Both are inter-exchangeable (symmetric measure).
# returns: a single real value between 0 and 1 indicating the proportion of elements
#		   put on the same part in both partitions.
############################################################################
process.percent.correct <- function(partition1, partition2)
{	# process the confusion matrix
	conf.matrix <- process.confusion.matrix(partition1,partition2)
	
	# init
	total <- 0
	norm <- sum(conf.matrix)
	conf.matrix <- conf.matrix / norm
	
	# match the communities
	while(nrow(conf.matrix)>0 && ncol(conf.matrix)>0)
	{	# get the position of the maximal value in the whole matrix
		mx <- which(conf.matrix == max(conf.matrix), arr.ind = TRUE)
		# add to score
		total <- total + conf.matrix[mx[1,1],mx[1,2]]
		# remove corresponding row & col
		conf.matrix <- conf.matrix[-mx[1,1],-mx[1,2],drop=FALSE]
	}
	
	# process final result
	return(total)
}



############################################################################
# A variant of the percent correct measure defined above, but this one is
# not symmetric: the first partition is compared to the second (which, here,
# is different from the opposite). 
# 
# For each part in partition #1, we look for the part in partition #2 with the
# largest intersection. These are considered to be matched. We discard them (both
# corresponding row and column are removed from the confusion matrix) and repeat 
# the same process until there is no more row to process.
#
# In Newman's stricter variant, which is activated by setting the no.merge option 
# to TRUE, if one part from partition #2 maximally intersects several parts from 
# partition #1, then *all* its elements are considered as misclassified.
#		Newman, M. E. J. 
#		Fast algorithm for detecting community structure in networks 
#		Physical Review E, 69:066133, 2004
#		http://pre.aps.org/abstract/PRE/v69/i6/e066133
# Also, note that in its seminal paper, Newman assesses the purity of the ground-
# truth partition relatively to the partition estimated by the considered community 
# detection algorithm.
#		Girvan, M. & Newman, M. E. J. 
#		Community structure in social and biological networks 
#		Proceedings of the National Academy of Sciences, 99:7821-7826, 2002
#		http://www.pnas.org/content/99/12/7821.abstract
# Most authors chose to do the opposite (i.e. estimated partition compared to 
# ground-truth), as explained in my own paper (see the README.md file for the 
# full reference).
#
# partition1: the first partition to consider, represented as an integer vector.
#			  Each value in the vector represents the id of a part. The parts 
#			  must be counted starting from one (not zero).
# partition2: the second partition to consider. Same representation as for
#			  the first one. Both are inter-exchangeable (symmetric measure).
# no.merge: if TRUE, applies Newman's stricter version of the measure.
# returns: a single real value between 0 and 1 indicating the purity of the first
#		   partition relatively to the second.
############################################################################
process.purity <- function(partition1, partition2, no.merge=FALSE)
{	# process the confusion matrix
	conf.matrix <- process.confusion.matrix(partition1,partition2)
	
	# init
	total <- 0
	norm <- sum(conf.matrix)
	conf.matrix <- conf.matrix / norm
	corresp <- rep(0,nrow(conf.matrix))
	
	# for each part in partition1, identify the corresponding parts in partition2
	# i.e. those with the largest intersection (there can be several)
	for(r in 1:nrow(conf.matrix))
		corresp[r] <- which.max(conf.matrix[r,])
	
	# identify and count the correctly classified nodes
	for(c in 1:ncol(conf.matrix))
	{	indices <- which(corresp==c)
		# possibly ignore merged communities
		if(length(indices)==1 || !no.merge)
			total <- total + sum(conf.matrix[indices,c])
	}
	
	return(total)
}



############################################################################
# Modification of the Purity measure described in my paper (cf. README.md for
# the full reference). It takes some topological information into account, in
# order to weight differently classification mistakes depending on the network
# position of the concerned nodes. The parameter weights must be a numerical
# vector countaining as many values as there are nodes in the processed network.
# It corresponds to the weights mentioned above.
#
# partition1: the first partition to consider, represented as an integer vector.
#			  Each value in the vector represents the id of a part. The parts 
#			  must be counted starting from one (not zero).
# partition2: the second partition to consider. Same representation as for
#			  the first one. Both are inter-exchangeable (symmetric measure).
# no.merge: if TRUE, applies Newman's stricter version of the measure.
# weights: numerical vector representing the topological weight associated
#		  to each node. So, its length must be the same as vectors partition1 
# 		  and partition2.
# returns: a single real value between 0 and 1 indicating the purity of the first
#		   partition relatively to the second.
############################################################################
process.topological.purity <- function(partition1, partition2, weights, no.merge=FALSE)
{	# process the confusion matrix
	conf.matrix <- process.weighted.confusion.matrix(partition1,partition2,weights)
	
	# init
	total <- 0
	norm <- sum(conf.matrix)
	conf.matrix <- conf.matrix / norm
	corresp <- rep(0,nrow(conf.matrix))
	
	# for each part in partition1, identify the corresponding parts in partition2
	# i.e. those with the largest intersection (there can be several)
	for(r in 1:nrow(conf.matrix))
		corresp[r] <- which.max(conf.matrix[r,])
	
	# identify and count the correctly classified nodes
	for(c in 1:ncol(conf.matrix))
	{	indices <- which(corresp==c)
		# possibly ignore merged communities
		if(length(indices)==1 || !no.merge)
			total <- total + sum(conf.matrix[indices,c])
	}
	
	return(total)
}
