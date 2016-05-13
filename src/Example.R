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




# set up some example partitions
part.ref <- c(1,1,1,1,1,1,2,2,2,2)		# ground truth
part.est <- list(
	AlgoA=c(1,2,1,1,1,1,2,2,2,2),		# estimation by algorithm A (paper)
	AlgoB=c(1,1,1,1,1,2,2,2,2,2),		# estimation by algorithm B (paper)
	Three=c(1,1,1,2,2,2,2,3,3,3),		# three parts
	Snglt=c(1,2,3,4,5,6,7,8,9,10),		# singletons
	Alli1=c(1,1,1,1,1,1,1,1,1,1),		# all in one
	Unbld=c(1,1,1,1,1,1,1,1,1,2)		# 9 vs. 1 (unbalanced)
)
# display these partitions
cat("Reference partition:");print(part.ref)
cat("Estimated partitions:\n")
for(e in 1:length(part.est))
{	cat("- ",names(part.est)[e],": ")
	print(part.est[[e]])
}



# create the example graph
g <- graph(edges=c(1,2,1,3,1,4,1,5,1,6,2,3,2,4,2,5,2,6,2,7,3,4,4,5,6,8,7,8,7,9,7,10,8,9,8,10,9,10), 
		directed=FALSE)
# plot graph and ground truth
plot(g, layout=layout.fruchterman.reingold, vertex.label=1:vcount(g)
	,vertex.color=part.ref+1
)



# process the topological weights
topo.weights <- list(
	NormalizedDegree=degree(g)/max(degree(g)),						# normalized degree
	Embeddedness=embeddedness(g,part.ref),							# embeddedness
	Product=degree(g)/max(degree(g))*embeddedness(g,part.ref)		# product of normalized degree and embeddedness
)
for(t in 1:length(topo.weights))
{	cat(names(topo.weights)[t],": ")
	print(topo.weights[[t]])
}



# display confusion matrices (for debug)
for(e in 1:length(part.est))
{	cat(names(part.est)[e],":Refrc\n",sep="")
	# confusion matrix
	cat(" - Confusion matrix - unweighted:\n")
	cm <- process.confusion.matrix(part.est[[e]],part.ref)
	print(cm)
	# weighted confusion matrices
	for(t in 1:length(topo.weights))
	{	cat(" - Confusion matrix - weights=",names(topo.weights)[[t]],":\n",sep="")
		cm <- process.weighted.confusion.matrix(part.est[[e]],part.ref,topo.weights[[t]])
		print(cm)
	}
	# pairwise comparison matrix
	cat(" - Pairwise comparison matrix - unweighted:\n")
	pm <- process.pairwise.matrix(part.est[[e]],part.ref)
	print(pm)
	# weighted pairwise comparison matrices
	for(t in 1:length(topo.weights))
	{	cat(" - Pairwise comparison matrix - weights=",names(topo.weights)[[t]],":\n",sep="")
		pm <- process.pairwise.matrix(part.est[[e]],part.ref,topo.weights[[t]])
		print(pm)
	}
}



# print result header
cat(sprintf("%-60s",""),"\tRefrc:Refrc\t",sep="")
for(e in 1:length(part.est))
	cat(names(part.est)[e],":Refrc\t",sep="")
cat("\n")



# list the traditional measures (no topological weight)
trad.measures <- list(
	Purity=function(p1,p2) process.purity(p1, p2, no.merge=FALSE),
	InversePurity=function(p1,p2) process.purity(p2, p1, no.merge=FALSE),
	NewmanPurity=function(p1,p2) process.purity(p2, p1, no.merge=TRUE),
	MeanPurity=function(p1,p2) harmonic.mean(process.purity(p1, p2, no.merge=FALSE), process.purity(p2, p1, no.merge=FALSE)),
	PercentCorrect=process.percent.correct,
	RandIndex=function(p1,p2) compare(p1, p2, method="rand"),
	AdjustedRandIndex=function(p1,p2) compare(p1, p2, method="adjusted.rand"),
	NormalizedMutualInformation=function(p1,p2) compare(p1, p2, method="nmi")
)
# process them
for(m in 1:length(trad.measures))
{	f <- trad.measures[[m]]
	val  <- f(part.ref, part.ref)
	cat(sprintf("%-60s\t%.9f",
		names(trad.measures)[m],val),sep="")
	for(e in 1:length(part.est))
	{	val <- f(part.est[[e]], part.ref)
		cat(sprintf("\t%.9f",val))
	}
	cat("\n")
}



# list the measures using a topological weight
topo.measures <- list(
	WeightedPurity=function(p1,p2,tm) process.topological.purity(p1, p2, weights=tm, no.merge=FALSE),
	WeightedInversePurity=function(p1,p2,tm) process.topological.purity(p2, p1, weights=tm, no.merge=FALSE),
	WeightedMeanPurity=function(p1,p2,tm) harmonic.mean(process.topological.purity(p1, p2, weights=tm, no.merge=FALSE), process.topological.purity(p2, p1, weights=tm, no.merge=FALSE)),
	WeightedRandIndex=function(p1,p2,tm) process.topological.rand.index(p1, p2, weights=tm),
	WeightedAdjustedRandIndex=function(p1,p2,tm) process.topological.adjusted.rand.index(p1, p2, weights=tm),
	WeightedNormalizedMutualInformation=function(p1,p2,tm) process.topological.NMI(p1, p2, weights=tm)
)
# process them
cat("------------------------\n")
for(m in 1:length(topo.measures))
{	f <- topo.measures[[m]]
	for(w in 1:length(topo.weights))
	{	val  <- f(part.ref, part.ref, topo.weights[[w]])
		cat(sprintf("%-60s\t%.9f",
			paste(names(topo.measures)[m],":",names(topo.weights)[w],sep=""),val),sep="")
		for(e in 1:length(part.est))
		{	val <- f(part.est[[e]], part.ref, topo.weights[[w]])
			cat(sprintf("\t%.9f",val))
		}
		cat("\n")
	}
}
