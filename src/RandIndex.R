############################################################################
# Functions used to process the Rand index and its variants.
# 
# Note: the Rand index is now implemented in igraph
# 		compare(partition1, partition2, method="rand")
# as well as its adjusted version:
# 		compare(partition1, partition2, method="adjusted.rand")
# It is also available in the library flexclust:
# 		randIndex(partition1,partition2,correct=FALSE)
# and the adjusted version:
# 		randIndex(partition1,partition2,correct=TRUE)
#
# Vincent Labatut 2012-16
############################################################################
# get the common functions
source("src/CommonFunctions.R")




############################################################################
# Modification of the Rand index allowing to use a different weight for each
# element, as described in my paper (cf. README.md for the full reference). 
#
# partition1: the first partition to consider, represented as an integer vector.
#			  Each value in the vector represents the id of a part. The parts 
#			  must be counted starting from one (not zero).
# partition2: the second partition to consider. Same representation than for
#			  the first one. Both are inter-exchangeable (symmetric measure).
# no.merge: if TRUE, applies Newman's stricter version of the measure.
# weights: numerical vector representing the topological weight associated
#		   to each node. So, its length must be the same than vectors partition1 
# 		   and partition2.
# returns: a single real value between 0 and 1 corresponding to the modified
#		   version of the Rand index.
############################################################################
process.topological.rand.index <- function(partition1, partition2, weights)
{	# init
	pwc.mat <- process.pairwise.matrix(partition1, partition2, weights)
	
	# process the measure
	result <- (pwc.mat[1,1]+pwc.mat[2,2]) / sum(pwc.mat)
	
	return(result)
}



############################################################################
# Processes the generalized binomial coefficient used in our version of the 
# adjusted Rand index calculation. The specified indices correspond to positions 
# in the original sequence of elements. The weights vector contains the weights
# for *all* nodes (not just those listed in the indices parameter).
#
# indices: ids of the elements considered in this calculation. These correspond
#		   to positions in the weights vector, and possiblu describe a subset
#		   of all the available elements.
# weights: individual weight associated with each element, including those not
#		   specificed in the indices parameter.
# returns: a real value corresponding to some kind of weighted version of the
#          binomial coefficient, in the case where we consider combinations of
#		   two elements in some set of n elements.
############################################################################
generalized.coefficient <- function(indices, weights)
{	# init
	size <- length(indices)
	result <- 0
	
	if(size>1)
	{	# consider each pair and their weight
		for(i in 1:(size-1))
		{	for(j in (i+1):size)
			{	idx1 <- indices[i]
				idx2 <- indices[j]
				result <- result + weights[idx1]*weights[idx2]
				#cat(weights[idx1]*weights[idx2],":",sep="")
			}
		}
	}
	
	return(result)
}



############################################################################
# Modification of the adjusted Rand index allowing to use a different weight 
# for each element, as described in my paper (cf. README.md for the full 
# reference). 
#
# partition1: the first partition to consider, represented as an integer vector.
#			  Each value in the vector represents the id of a part. The parts 
#			  must be counted starting from one (not zero).
# partition2: the second partition to consider. Same representation than for
#			  the first one. Both are inter-exchangeable (symmetric measure).
# no.merge: if TRUE, applies Newman's stricter version of the measure.
# weights: numerical vector representing the topological weight associated
#		   to each node. So, its length must be the same than vectors partition1 
# 		   and partition2.
# returns: a single real value between 0 and 1 corresponding to the modified
#		   version of the adjusted Rand index.
############################################################################
process.topological.adjusted.rand.index<-function(partition1, partition2, weights)
{	# whole graph
	ws <- generalized.coefficient(1:length(partition1),weights)
	
	# parts 1
	s1 <- max(partition1)
	w1 <- 0
	for(i in 1:s1)
	{	indices1 <- which(partition1==i)
		w1 <- w1 + generalized.coefficient(indices1,weights)
	}
	
	# parts 2
	s2 <- max(partition2)
	w2 <- 0
	for(j in 1:s2)
	{	indices2 <- which(partition2==j)
		w2 <- w2 + generalized.coefficient(indices2,weights)
	}
	
	# intersections
	w12 <- 0
	for(i in 1:s1)
	{	indices1 <- which(partition1==i)
		for(j in 1:s2)
		{	indices2 <- which(partition2==j)
			indices12 <- intersect(indices1,indices2)
			w12 <- w12 + generalized.coefficient(indices12,weights)
		}
	}
	
	# process measure (cf. the formula in the paper)
	result <- (w12 - w1*w2/ws) / (1/2*(w1+w2) - w1*w2/ws)
	return(result)
}
