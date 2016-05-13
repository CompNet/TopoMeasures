############################################################################
# Functions used to process the normalized mutual information and its variants.
# 
# Note: the NMI is now implemented in igraph
# 		compare(partition1, partition2, method="nmi")
#
# Vincent Labatut 2012-16
############################################################################
# get the common functions
source("src/CommonFunctions.R")




############################################################################
# There are several different variant of the mutual information called "normalized",
# which differ in the way this normalization is performed, see for instance:
#		https://en.wikipedia.org/wiki/Mutual_information#Normalized_variants
# Here, we use the version introduced by Danon et al. to the field of community
# detection, i.e. 2I(X,Y)/(H(X)+H(Y)), as described in:
#		Danon, L.; Díaz-Guilera, A.; Duch, J. & Arenas, A. 
#		Comparing Community Structure Identification 
#		Journal of Statistical Mechanics, 9:P09008, 2005
# 
# partition1: the first partition to consider, represented as an integer vector.
#			  Each value in the vector represents the id of a part. The parts 
#			  must be counted starting from one (not zero).
# partition2: the second partition to consider. Same representation than for
#			  the first one. Both are inter-exchangeable (symmetric measure).
# weights: numerical vector representing the topological weight associated
#		   to each node. So, its length must be the same than vectors partition1 
# 		   and partition2.
# returns: a single real value between 0 and 1 corresponding to the NMI version
#		   used in the community detection field.
############################################################################
process.topological.NMI <- function(partition1, partition2, weights)
{	# process the confusion matrix
	conf.matrix <- process.weighted.confusion.matrix(partition1,partition2,weights)
	
	# init
	norm <- sum(conf.matrix)
	conf.matrix <- conf.matrix / norm
	
	# process all sums once and for all
	s.row <- rowSums(conf.matrix)
	s.col <- colSums(conf.matrix)
	
	# process normalized mutual information measure
	sum1 <- 0
	for(i in 1:nrow(conf.matrix))
	{	for(j in 1:ncol(conf.matrix))
		{	temp <- conf.matrix[i,j] * log2((conf.matrix[i,j]) / (s.row[i]*s.col[j]))
			if(!is.na(temp))
				sum1 <- sum1 + temp
		}
	}
	
	sum2 <- 0
	for(i in 1:nrow(conf.matrix))
	{	temp <- s.row[i]*log2(s.row[i])
		if(!is.na(temp))
			sum2 <- sum2 + temp
	}
	
	sum3 <- 0
	for(j in 1:ncol(conf.matrix))
	{	temp <- s.col[j]*log2(s.col[j])
		if(!is.na(temp))
			sum3 <- sum3 + temp
	}
	
	# combine all the parts...
	result <- (-2)*sum1 / (sum2+sum3)
	return(result)
}
