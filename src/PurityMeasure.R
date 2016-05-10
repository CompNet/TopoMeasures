#####################################################################
# Functions used to process the purity measure and its variants.
#####################################################################






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
