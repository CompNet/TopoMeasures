############################################################################
# This script shows how to use the proposed measures, and compares their 
# results with that of traditional measures, on a few examples including
# the graph from the paper (cf. README.md for the full reference).
#
# Vincent Labatut 2012-16
############################################################################

#setwd("C:/Eclipse/workspaces/Networks/TopoMeasures")
#setwd("~/eclipse/workspaces/Networks/TopoMeasures")
#source("src/Example.R")

# get the measure scripts
source("src/NormalizedMutualInformation.R")
source("src/PurityMeasure.R")
source("src/RandIndex.R")



# create the example graph
g <- graph(edges=c(1,2,1,3,1,4,1,5,1,6,2,3,2,4,2,5,2,6,2,7,3,4,4,5,6,8,7,8,7,9,7,10,8,9,8,10,9,10), 
		directed=FALSE)
# set up some example partition
ref <- c(1,1,1,1,1,1,2,2,2,2)		# ground truth
est1 <- c(1,1,1,1,1,2,2,2,2,2)		# estimation by algorithm #1
est2 <- c(1,2,1,1,1,1,2,2,2,2)		# estimation by algorithm #2
# plot graph and ground truth
plot(g, layout=layout.fruchterman.reingold, vertex.label=1:vcount(g)
	,vertex.color=ref+1
)


# process the nodal measures later used as weights
deg <- degree(g)
emb <- embeddedness(g)


# list the measures we want to process for this graph and partitions
measures <- list(	
	# traditional measures (no topological weight)
	Purity=function(p1,p2) process.purity(p1, p2, no.merge=FALSE),
	InversePurity=function(p1,p2) process.purity(p2, p1, no.merge=FALSE),
	NewmanPurity=function(p1,p2) process.purity(p2, p1, no.merge=TRUE),
	MeanPurity=function(p1,p2) harmonic.mean(process.purity(p1, p2, no.merge=FALSE), process.purity(p2, p1, no.merge=FALSE)),
	PercentCorrect=process.percent.correct,
	RandIndex=function(p1,p2) compare(p1, p2, method="rand"),
	AdjustedRandIndex=function(p1,p2) compare(p1, p2, method="adjusted.rand"),
	NormalizedMutualInformation=function(p1,p2) compare(p1, p2, method="nmi"),
	
	# using degree as the topological weight
	WeightedDegreePurity=function(p1,p2) process.topological.purity(p1, p2, topo.measure=deg, no.merge=FALSE,),
	WeightedDegreeInversePurity=function(p1,p2) process.topological.purity(p2, p1, no.merge=FALSE),
	WeightedDegreeMeanPurity=function(p1,p2) harmonic.mean(process.topological.purity(p1, p2, topo.measure=deg, no.merge=FALSE), process.topological.purity(p2, p1, topo.measure=deg, no.merge=FALSE)),
	WeightedDegreeRandIndex=function(p1,p2) process.topological.rand.index(p1, p2, topo.measure=deg),
	WeightedDegreeAdjustedRandIndex=function(p1,p2) process.topological.adjusted.rand.index(p1, p2, topo.measure=deg),
	WeightedDegreeNormalizedMutualInformation=function(p1,p2) process.topological.NMI(p1, p2, topo.measure=deg),

	# using embeddedness as the topological weight
	WeightedDegreePurity=function(p1,p2) process.topological.purity(p1, p2, topo.measure=emb, no.merge=FALSE,),
	WeightedDegreeInversePurity=function(p1,p2) process.topological.purity(p2, p1, no.merge=FALSE),
	WeightedDegreeMeanPurity=function(p1,p2) harmonic.mean(process.topological.purity(p1, p2, topo.measure=emb, no.merge=FALSE), process.topological.purity(p2, p1, topo.measure=emb, no.merge=FALSE)),
	WeightedDegreeRandIndex=function(p1,p2) process.topological.rand.index(p1, p2, topo.measure=emb),
	WeightedDegreeAdjustedRandIndex=function(p1,p2) process.topological.adjusted.rand.index(p1, p2, topo.measure=emb),
	WeightedDegreeNormalizedMutualInformation=function(p1,p2) process.topological.NMI(p1, p2, topo.measure=emb),
	
	# using degree*embeddedness as the topological weight
	WeightedDegreePurity=function(p1,p2) process.topological.purity(p1, p2, topo.measure=deg*emb, no.merge=FALSE,),
	WeightedDegreeInversePurity=function(p1,p2) process.topological.purity(p2, p1, no.merge=FALSE),
	WeightedDegreeMeanPurity=function(p1,p2) harmonic.mean(process.topological.purity(p1, p2, topo.measure=deg*emb, no.merge=FALSE), process.topological.purity(p2, p1, topo.measure=deg*emb, no.merge=FALSE)),
	WeightedDegreeRandIndex=function(p1,p2) process.topological.rand.index(p1, p2, topo.measure=deg*emb),
	WeightedDegreeAdjustedRandIndex=function(p1,p2) process.topological.adjusted.rand.index(p1, p2, topo.measure=deg*emb),
	WeightedDegreeNormalizedMutualInformation=function(p1,p2) process.topological.NMI(p1, p2, topo.measure=deg*emb)
)


# process each measure
for(m in 1:length(measures))
{	f <- mesures[[m]]
	sr  <- f(ref, ref)
	se1 <- f(se1, ref)
	se2 <- f(se2, ref)
	cat(names(measures)[m],": ","ref vs. ref:",sr,"\test1 vs. ref:",se1,"\test2 vs ref:",se2,"\n",sep="")
}
