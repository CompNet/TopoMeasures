######################################################################
# these functions are related to the processing of performance measures
# (in the context of community detection)
# source("CommunityDetection/performances/processing.R")
######################################################################

#library("mclust") #TODO is this really needed ?




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
