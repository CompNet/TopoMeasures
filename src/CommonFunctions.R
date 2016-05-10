#####################################################################
# Functions used in the processing of more than one measure.
#####################################################################
# main lib, used to represent and process graphs
library("igraph")




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
	
	# process the matrix manually
#	conf.matrix <- matrix(0,nrow=num.row,ncol=num.col)
#	for(i in 1:length(part1))
#	{	row <- part1[i]
#		col <- part2[i]
#		conf.matrix[row,col] <- conf.matrix[row,col] + 1
#	}
	
	
	# more efficient: use the standard 'table' function
	# note: some rows/cols can be missing if they don't appear in at least one part
	conf.matrix <- as.matrix(table(part1,part2)) 
	
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



############################################################################
# Processes the embeddedness as defined by Lancichinetti et al. in:
# 		Lancichinetti, A.; Kivelä, M.; Saramäki, J. & Fortunato, S. 
#		Characterizing the Community Structure of Complex Networks 
#		PLoS ONE, 2010, 5, e11976
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
		#cat("com1=",com1,"\n")
		
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
