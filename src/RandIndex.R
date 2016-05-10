#####################################################################
# Functions used to process the Rand index and its variants.
#####################################################################
# fast implementation of the Rand index and its adjusted version
library("flexclust") 

# get the common functions
source("src/CommonFunctions.R")



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
# the rand index is nothing but the Cohen's Kappa statistic 
# processed on the pairwise matrix
# TODO DONE
process.adjusted.rand.index<-function(part1, part2)
{	#conf.matrix <- process.confusion.matrix(part1,part2);
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
