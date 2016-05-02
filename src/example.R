#setwd("c:/eclipse/workspaces/Networks")
setwd("~/eclipse/workspaces/Networks")

#source("c:/Documents and Settings/Vincent/Mes documents/Travail/Ecrits/Networks Analysis/Papiers/2012.3.MARAMI12 extended/_1.scripts/figure1.R")
#source("c:/Eclipse/workspaces/Networks/CommunityDetection/performances/processing.R")
source("CommunityDetection/performances/processing.R")

#source("Process/example.R")

library(igraph)


#g <- graph(edges=c(
#	1,2,1,3,1,4,1,5,1,6,2,3,2,4,2,5,2,6,2,7,3,4,4,5,6,8,7,8,7,9,7,10,8,9,8,10,9,10
#)-1, directed=FALSE) # for versions of igraph where numbering starts from 0
g <- graph(edges=c(
	1,2,1,3,1,4,1,5,1,6,2,3,2,4,2,5,2,6,2,7,3,4,4,5,6,8,7,8,7,9,7,10,8,9,8,10,9,10
), directed=FALSE)

ref <- c(1,1,1,1,1,1,2,2,2,2)
est1 <- c(1,1,1,1,1,2,2,2,2,2)
est2 <- c(1,2,1,1,1,1,2,2,2,2)

plot(g, layout=layout.fruchterman.reingold, vertex.label=1:vcount(g)
	,vertex.color=ref+1
)

measures <- c("MODULARITY",

		"NORM_MUT_INF",

		"RAND_IDX",
		"ADJ_RAND_IDX",

		"PRCNT_CRCT",
		"PURITY_REF",
		"PURITY_EST",
		"PURITY_BOTH",

		"TOPO_PURITY_REF_DEGREE",
		"TOPO_PURITY_REF_EMBEDDEDNESS",
		"TOPO_PURITY_REF_BOTH",
		"TOPO_PURITY_EST_DEGREE",
		"TOPO_PURITY_EST_EMBEDDEDNESS",
		"TOPO_PURITY_EST_BOTH",
		"TOPO_PURITY_BOTH_DEGREE",
		"TOPO_PURITY_BOTH_EMBEDDEDNESS",
		"TOPO_PURITY_BOTH_BOTH",

		"TOPO_ADJ_RAND_IDX_DEGREE",
		"TOPO_ADJ_RAND_IDX_EMBEDDEDNESS",
		"TOPO_ADJ_RAND_IDX_BOTH",
		"TOPO_RAND_IDX_DEGREE",
		"TOPO_RAND_IDX_EMBEDDEDNESS",
		"TOPO_RAND_IDX_BOTH",

		"TOPO_NORM_MUT_INF_DEGREE",
		"TOPO_NORM_MUT_INF_EMBEDDEDNESS",
		"TOPO_NORM_MUT_INF_BOTH"
)
for(measure in measures)
{	x <- process_performance_metric(g=g, partition=ref, reference=ref, perfname=measure)
	y <- process_performance_metric(g=g, partition=est1, reference=ref, perfname=measure)
	z <- process_performance_metric(g=g, partition=est2, reference=ref, perfname=measure)
	cat(measure,": ","R=",x," E1=",y," E2=",z,"\n",sep="")
}

# figure 2
#g2 <- graph.star(n=8,mode="undirected")
#plot(g2,layout=layout.star(g2))
