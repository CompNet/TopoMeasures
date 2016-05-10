#####################################################################
# Functions used to process the normalized mutual information and its variants.
#####################################################################



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
