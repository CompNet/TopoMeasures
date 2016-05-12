############################################################################
# Functions used in the processing of more than one measure.
############################################################################
# main lib, used to represent and process graphs
library("igraph")




############################################################################
# Processes the confusion matrix for the specified partitions, i.e. the matrix 
# whose cell (i,j) contains the number of elements belonging to part #i in 
# partition #1 and to part #j in partition #2.
# 
# partition1: the first partition to consider, represented as an integer vector.
#			  Each value in the vector represents the id of a part. The parts 
#			  must be counted starting from one (not zero).
# partition2: the second partition to consider. Same representation than for
#			  the first one. Both are inter-exchangeable (symmetric measure).
# returns: the confusion matrix, a.k.a. contingency table
############################################################################
process.confusion.matrix <- function(part1, part2)
{	p1 <- factor(part1, levels=1:max(part1))
	p2 <- factor(part2, levels=1:max(part2))
	result <- as.matrix(table(p1,p2))
	return(result)
}



############################################################################
# Processes the pairwise comparison matrix for the specified partitions, i.e. 
# the 2x2 matrix containing the following values:
#	- (1,1): number of pairs of elements belonging to the same part in both partitions.
#	- (2,2): number of pairs of elements belonging to different parts in both partitions.
#	- (1,2),(2,1): number of pairs of elements belonging to different parts in one partition, 
#				   and to the same part in the other.
#
# partition1: the first partition to consider, represented as an integer vector.
#			  Each value in the vector represents the id of a part. The parts 
#			  must be counted starting from one (not zero).
# partition2: the second partition to consider. Same representation than for
#			  the first one. Both are inter-exchangeable (symmetric measure).
############################################################################
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
					#cat("same-same\n")				
					result[1,1] <- result[1,1] + 1
				}
				else
				{	# the same for est., different for ref.
					#cat("same-diff\n")				
					result[1,2] <- result[1,2] + 1
				}
			}
			else
			{	if(part2[i]==part2[j])
				{	# different for est., the same for ref.
					#cat("diff-same\n")				
					result[2,1] <- result[2,1] + 1
				}
				else
				{	# different for both
					#cat("diff-diff\n")				
					result[2,2] <- result[2,2] + 1
				}
			}
		}
	}
	
	return(result)
}



############################################################################
# Processes the embeddedness as defined by Lancichinetti et al. in:
# 		Lancichinetti, A.; Kivelä, M.; Saramäki, J. & Fortunato, S. 
#		Characterizing the Community Structure of Complex Networks 
#		PLoS ONE, 5:e11976, 2010
#		http://www.plosone.org/article/info:doi/10.1371/journal.pone.0011976
# It is basically defined as k_{in}/k, i.e. the ratio of the internal to the
# total degree of the considered node. The internal degree is its number of
# connections inside its own community. The total degree is its total number
# of connections (ignoring the neighbors' communities).
# 
# g: graph to process.
# partition: integer vector of vcount(g) values, representing a partition of
#			 the graph node set. Each value corresponds to the id of the part 
#			 containing the corresponding node. These id must start from one,
#			 not zero.
# returns: a vector of vcount(g) numerical values, each one corresponding to
#		   the embeddedness of a node.
############################################################################
embeddedness <- function(g, partition)
{	result <- rep(NA,vcount(g))
	for(idx1 in 1:vcount(g))
	{	# init
		com1 <- partition[idx1]
		int.deg <- 0
		ext.deg <- 0
		
		# get the direct neighbors
		neighs <- neighborhood(graph=g,order=1,nodes=idx1,mode="all")[[1]]
		neighs <- neighs[neighs!=idx1]
		#print(neighs)
		#print(partition[neighs])
		
		# check if they are in the same community
		for(idx2 in neighs)
		{	com2 <- partition[idx2]
			if(com1==com2)
				int.deg <- int.deg + 1
			else
				ext.deg <- ext.deg + 1
		}
		
		# process the final value
		result[idx1] <- int.deg/(int.deg+ext.deg)
		#cat("internal=",int.deg," external=",ext.deg,"\n")		
		#cat("degree=",degree(g)[i]," estimation=",int.deg+ext.deg,"\n")		
	}
	
	return(result)
}
