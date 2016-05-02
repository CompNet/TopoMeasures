############################################################################
# setwd("~/eclipse/workspaces/Networks")
# source("Process/gunce.comdet.R")
############################################################################
library("igraph")

source("CommonTools/networks/cleaning.R")
source("CommonTools/networks/converting.R")
source("CommonTools/networks/linegraph.R");
source("CommonTools/networks/plotting.R");
source("CommonTools/parameters/correlating.R")
source("CommonTools/parameters/distribution_plotting.R")
source("CommonTools/parameters/overall_single_plotting.R")
source("CommonTools/parameters/overall_multiple_plotting.R")
source("CommonTools/parameters/files.R")
source("CommonTools/properties/aggregating.R")
source("CommonTools/properties/general.R")
#source("CommonTools/properties/processing.R")
source("CommonTools/test/generating.R")

source("CommunityDetection/agreement/processing.R")
source("CommunityDetection/algorithms/general.R");
source("CommunityDetection/partitions/aggregating.R")
source("CommunityDetection/partitions/distributions.R")
source("CommunityDetection/performances/averaging.R")
source("CommunityDetection/performances/reprocessing.R")
source("CommunityDetection/processing/byalgorithm.R")
source("CommunityDetection/processing/byiteration.R")
source("CommunityDetection/plots/comparisons.R")
source("CommunityDetection/plots/distributions.R")
source("CommunityDetection/plots/hierarchical.R")
source("CommunityDetection/plots/parameters.R")




#################### init variables ################### 
cat(":: init variables\n")
#folder <- "/home/bkantarci/Desktop/NewmanTransitive/n1000"							# TODO unix only
folder <- "/var/data/deneme2"							# TODO unix only
#folder <- "/media/DABE5366BE533A69/Backup/var/data/deneme2"							# TODO unix only
#folder <- "c:/data"								# TODO windows only
format <- "PNG"									# format for the generated plots
layout <- "FR"									# layout algorithm for the generated network plots
directed <- FALSE								# nature of the links in the generated networks
weighted <- FALSE								# presence of weights in the generated networks
type <- "pajek"									# type of file used to store the generated networks
prefix <- NA									# file prefix
verbose <- 10									# verbosity level
iterations <- 0:4  								# iterations
parameters_names <- c("N","AVGRDEGREE","GAMMA","BETA","MU")		# parameters internal names
parameters_texts <- c("n","avrgdegree","gamma","beta","mixing")	# parameters file names
#  parameters_values <- list(N=c(7500),AVGRDEGREE=c(6),GAMMA=c(3),BETA=c(2), MU=c(0))
  parameters_values <- list(N=c(25000),AVGRDEGREE=c(4),GAMMA=c(3),BETA=c(2), MU=c(0.2))
#  parameters_values <- list(N=c(250000),AVGRDEGREE=c(2),GAMMA=c(3),BETA=c(2), MU=c(0.2))
parameters_ranges <- list(N=c(7500,250000),AVGRDEGREE=c(2,6),GAMMA=c(3,3),BETA=c(2,2),MU=c(0,1))
parameters_data <- list(names=parameters_names,values=parameters_values,texts=parameters_texts,ranges=parameters_ranges)
algorithms <- c(
		"COPRA",
		"FASTGREEDY",
		"INFOMAP",
		"INFOMOD",
#		"LABELPROP",
		"LOUVAIN",
		"MARKOVCLUSTER",
		"OSLOM",
#		"RADICCHI",
		"WALKTRAP"
	)
performances_names <- c(						# names of the performances to be processed
#		"MODULARITY"		
#		"NORM_MUT_INF",
#		"PRCNT_CRCT",
#		"RAND_IDX",
		"ADJ_RAND_IDX",
#		"PURITY_REF",
#		"PURITY_EST",
#		"PURITY_BOTH",
#		"TOPO_PURITY_REF_DEGREE",
#		"TOPO_PURITY_REF_EMBEDDEDNESS",
#		"TOPO_PURITY_REF_BOTH",
#		"TOPO_PURITY_EST_DEGREE",
#		"TOPO_PURITY_EST_EMBEDDEDNESS",
#		"TOPO_PURITY_EST_BOTH",
#		"TOPO_PURITY_BOTH_DEGREE",
#		"TOPO_PURITY_BOTH_EMBEDDEDNESS",
#		"TOPO_PURITY_BOTH_BOTH",
#		"TOPO_PURITY_BOTH_CLOSENESS",
#		"TOPO_RAND_IDX_DEGREE",
#		"TOPO_RAND_IDX_EMBEDDEDNESS",
#		"TOPO_RAND_IDX_BOTH",
		"TOPO_ADJ_RAND_IDX_DEGREE",
		"TOPO_ADJ_RAND_IDX_EMBEDDEDNESS",
		"TOPO_ADJ_RAND_IDX_BOTH"
#		"TOPO_NORM_MUT_INF_DEGREE",
#		"TOPO_NORM_MUT_INF_EMBEDDEDNESS",
#		"TOPO_NORM_MUT_INF_BOTH"
)
#performances_names <- c("NORM_MUT_INF")
#######################################################



######### apply community detection algorithms ########
cat(":: apply community detection algorithms\n")
##	apply separately each algorithm in the specified list
#for(algoname in algorithms)
#apply_algo_to_parameters(folder,prefix,parameters_data,iterations,algoname,directed,weighted,performances_names,force=TRUE,show_console=FALSE,save_hierarchy=TRUE,verbose)
##	alternatively: apply all specified algorithms iteration by iteration 
#apply_algos_to_parameters(folder,prefix,parameters_data,iterations,algonames=c(algorithms),directed,weighted,performances_names,force=TRUE,show_console=TRUE,save_hierarchy=FALSE,verbose)
#######################################################



############### reprocess performances ################
#cat(":: reprocess performances\n")
##  Reprocess all performance measures from some previously estimated partitions.
##  Usefull if some bug is detected in the way performance is processed, 
##  or if some new performance measure is implemented
reprocess_performance_parameters(folder,prefix,parameters_data,iterations,algonames=algorithms,performances_names,delete=FALSE,save_hierarchy=FALSE,remove_singles=TRUE,verbose)
#######################################################



################ average performances #################
#cat(":: average performances\n")
##	Average the performances either by using everything in the main folder or the specified parameters
average_performance_parameters(folder,prefix,parameters_data,iterations,algonames=algorithms,delete=TRUE,verbose)
#######################################################



################### plot networks #####################
#cat(":: plot networks\n")
##  Generate the network plots.
##  Use all the estimated partitions (i.e. communities) to color the graphs.
##  In other term, one plot is generated for the reference partition, and one plot for each algorithm.
##  Note the exact same layout is used to ease visual comparisons.
#plot_networks_from_parameters(folder,prefix,algonames=algorithms,directed,parameters_data,iterations,filetype=type,layout,layout_level=0,format,force=TRUE,verbose)
#######################################################



################# plot performances ###################
#cat(":: plot performances\n")
##  Generate the performances plots,
##	i.e. lines plots comparing the performances depending on the parameters and algorithms
#plot_performances_parameters(folder,prefix,parameters_data,algonames=algorithms,performances_names,format,verbose)
##	plot the hierarchical results (only for hierarchical algos, of course)
## 	(if they have been previously processed)
#plot_hierarchical_from_parameters(folder,prefix,algonames=algorithms,performances_names,parameters_data,iterations,levels=FALSE,format,force=TRUE,verbose)
#######################################################



############ community size distribution ##############
#cat(":: community sizes distribution\n")
##	process the community sizes distribution
#process_community_distribution_parameters(folder,prefix,parameters_data,iterations,algonames=algorithms,delete=TRUE,verbose)
## 	sum them over iterations
#aggregate_community_distributions_parameters(folder,prefix,parameters_data,iterations,algonames=algorithms,delete=TRUE,verbose)
##	generate the corresponding plots
#distribution_plot_communities_parameters(folder,prefix,parameters_data,algonames=algorithms,format,verbose)
#######################################################
	


################ partition agreement ##################
#cat(":: process partition agreement\n")
##	process the agreement between the partitions estimated by the various algorithms
#partition_agreement_fixed_network_parameters(folder,prefix,parameters_data,iterations,algonames=algorithms,performances_names=c("NORM_MUT_INF","RAND_IDX","ADJ_RAND_IDX"),consider_hierarchy=FALSE,force=TRUE,verbose)
#######################################################



################ plot performances ##################
cat(":: plot performances\n")
#t <- read.table("/var/data/deneme2/n7500/avrgdegree6/gamma3/beta2/mixing0/data.performances")
t <- read.table("/var/data/deneme2/n25000/avrgdegree4/gamma3/beta2/mixing0.2/data.performances")
m <- as.matrix(t)
indices <- c()
for(name in c(#"PRCNT_CRCT_MEAN","PRCNT_CRCT_STDEV"
#				,"NORM_MUT_INF_MEAN","NORM_MUT_INF_STDEV"
#				,"PRCNT_CRCT_MEAN","PRCNT_CRCT_STDEV"
#				,"RAND_IDX_MEAN","RAND_IDX_STDEV"
				"ADJ_RAND_IDX_MEAN","ADJ_RAND_IDX_STDEV"
#				,"PURITY_REF_MEAN","PURITY_REF_STDEV"
#				,"PURITY_EST_MEAN","PURITY_EST_STDEV"
#				,"PURITY_BOTH_MEAN","PURITY_BOTH_STDEV"
#				,"TOPO_PURITY_REF_DEGREE_MEAN","TOPO_PURITY_REF_DEGREE_STDEV"
#				,"TOPO_PURITY_REF_EMBEDDEDNESS_MEAN","TOPO_PURITY_REF_EMBEDDEDNESS_STDEV"
#				,"TOPO_PURITY_REF_BOTH_MEAN","TOPO_PURITY_REF_BOTH_STDEV"
#				,"TOPO_PURITY_EST_DEGREE_MEAN","TOPO_PURITY_EST_DEGREE_STDEV"
#				,"TOPO_PURITY_EST_EMBEDDEDNESS_MEAN","TOPO_PURITY_EST_EMBEDDEDNESS_STDEV"
#				,"TOPO_PURITY_EST_BOTH_MEAN","TOPO_PURITY_EST_BOTH_STDEV"
#				,"TOPO_PURITY_BOTH_DEGREE_MEAN","TOPO_PURITY_BOTH_DEGREE_STDEV"
#				,"TOPO_PURITY_BOTH_EMBEDDEDNESS_MEAN","TOPO_PURITY_BOTH_EMBEDDEDNESS_STDEV"
#				,"TOPO_PURITY_BOTH_BOTH_MEAN","TOPO_PURITY_BOTH_BOTH_STDEV"
#				,"TOPO_PURITY_BOTH_CLOSENESS_MEAN","TOPO_PURITY_BOTH_CLOSENESS_STDEV"
#				,"TOPO_RAND_IDX_DEGREE_MEAN","TOPO_RAND_IDX_DEGREE_STDEV"
#				,"TOPO_RAND_IDX_EMBEDDEDNESS_MEAN","TOPO_RAND_IDX_EMBEDDEDNESS_STDEV"
#				,"TOPO_RAND_IDX_BOTH_MEAN","TOPO_RAND_IDX_BOTH_STDEV"
#				,"TOPO_ADJ_RAND_IDX_DEGREE_MEAN","TOPO_ADJ_RAND_IDX_DEGREE_STDEV"
#				,"TOPO_ADJ_RAND_IDX_EMBEDDEDNESS_MEAN","TOPO_ADJ_RAND_IDX_EMBEDDEDNESS_STDEV"
				,"TOPO_ADJ_RAND_IDX_BOTH_MEAN","TOPO_ADJ_RAND_IDX_BOTH_STDEV"
#				,"TOPO_NORM_MUT_INF_DEGREE_MEAN","TOPO_NORM_MUT_INF_DEGREE_STDEV"
#				,"TOPO_NORM_MUT_INF_EMBEDDEDNESS_MEAN","TOPO_NORM_MUT_INF_EMBEDDEDNESS_STDEV"
#				,"TOPO_NORM_MUT_INF_BOTH_MEAN","TOPO_NORM_MUT_INF_BOTH_STDEV"
			))
	indices <- c(indices,which(colnames(m)==name))
m <- m[,indices]
pdf(file="/home/vlabatut/figure.arj.pdf",bg="white",width=11)
x <- barplot(m[,seq(1,ncol(m)-1,2)],beside=TRUE, ylim=c(0,1),
	col=rainbow(nrow(m),v=0.8),
	legend.text=TRUE
)
for(i in 1:nrow(m))
	plotCI(x=x[i,], y=m[i,seq(1,ncol(m)-1,2)],
		uiw=m[i,seq(2,ncol(m),2)],
		pch=NA,err="y",add=TRUE)
dev.off()

################ plot details ##################
#reference <- load_partition(filename="/var/data/deneme2/n25000/avrgdegree4/gamma3/beta2/mixing0.2/it0/reference.clu")
#estimation <- load_partition(filename="/var/data/deneme2/n25000/avrgdegree4/gamma3/beta2/mixing0.2/it0/markovcluster.clu")
#g <- read.graph("/var/data/deneme2/n25000/avrgdegree4/gamma3/beta2/mixing0.2/it0/pajek.net",format="pajek")
##topo_measure <- degree(g)/max(degree(g))
##topo_measure <- embeddedness(g,reference)
##topo_measure <- degree(g)/max(degree(g))*embeddedness(g,reference)
##topo_measure <- sqrt(degree(g)/max(degree(g))*embeddedness(g,reference))
##topo_measure <- (degree(g)/max(degree(g))+embeddedness(g,reference))/2
##topo_measure <- 2*degree(g)/max(degree(g))*embeddedness(g,reference)/(degree(g)/max(degree(g))+embeddedness(g,reference))
##topo_measure <- betweenness(graph=g,directed=FALSE)
#topo_measure <- closeness(graph=g)
#distribution <- process_topological_purity(reference=reference, estimation=estimation, inverted=FALSE, topo_measure=topo_measure, details=TRUE)
##hist(topo_measure/sum(topo_measure))
#h <- hist(topo_measure)
#dev.new()
#h0 <- hist(topo_measure[distribution==0],breaks=h$breaks)
#dev.new()
#h1 <- hist(topo_measure[distribution==1],breaks=h$breaks)
##indices <- 2:length(h$counts)
##barplot(rbind(h1$counts[indices],h0$counts[indices]),col=c("blue","red"),names.arg=h$breaks[indices]
##	,log="y"
##)

################ statistical tests by algos ##################
## load individual results
#measure <- "TOPO_RAND_IDX_BOTH" # PURITY_BOTH TOPO_PURITY_BOTH_BOTH NORM_MUT_INF TOPO_NORM_MUT_INF_BOTH RAND_IDX TOPO_RAND_IDX_BOTH ADJ_RAND_IDX TOPO_ADJ_RAND_IDX_BOTH
#cat("measure:",measure,"\n")
#data <- matrix(nrow=length(algorithms),ncol=5)
#rownames(data) <- algorithms
#for(it in 0:4)
#{	#t <- read.table(paste("/var/data/deneme2/n7500/avrgdegree6/gamma3/beta2/mixing0/it",it,"/data.performances",sep=""))
#	t <- read.table(paste("/var/data/deneme2/n25000/avrgdegree4/gamma3/beta2/mixing0.2/it",it,"/data.performances",sep=""))
#	for(a in algorithms)
#	{	data[a,it+1] <- t[a,measure]
#	}
#}
## create factors
#data <- as.vector(t(data))
#groups <- factor(rep(algorithms,each=5))
## test for homogeneity of variances (homoskedasticity): p>0.05 means ok
#print(bartlett.test(data,groups))	# primary test
#print(fligner.test(data,groups))	# alternative test
## process anova: p>0.05 means can't say one mean is significantly different
#print(anova(lm(formula= data~groups)))
## post-hoc tests to identify which means are differents
#print(pairwise.t.test(data,groups,p.adj="none"))	# t-test without adjustment
#print(pairwise.t.test(data,groups,p.adj="bonf"))	# bonferroni adjustment
#print(pairwise.t.test(data,groups,p.adj="holm"))	# holm adjustment
#av <- aov(formula= data~groups)
#print(av)
#print(TukeyHSD(av))


################ statistical tests by measure ##################
# load individual results
#classic_measure <- "RAND_IDX" # PURITY_BOTH NORM_MUT_INF RAND_IDX
#topo_measure <- "TOPO_RAND_IDX_BOTH" # TOPO_PURITY_BOTH_BOTH TOPO_NORM_MUT_INF_BOTH TOPO_RAND_IDX_BOTH
#cat("classic measure:",classic_measure,"topological measure:",topo_measure,"\n")
#classic_data <- matrix(nrow=length(algorithms),ncol=5)
#rownames(classic_data) <- algorithms
#topo_data <- matrix(nrow=length(algorithms),ncol=5)
#rownames(topo_data) <- algorithms
#for(it in 0:4)
#{	#t <- read.table(paste("/var/data/deneme2/n7500/avrgdegree6/gamma3/beta2/mixing0/it",it,"/data.performances",sep=""))
#	t <- read.table(paste("/var/data/deneme2/n25000/avrgdegree4/gamma3/beta2/mixing0.2/it",it,"/data.performances",sep=""))
#	for(a in algorithms)
#	{	classic_data[a,it+1] <- t[a,classic_measure]
#		topo_data[a,it+1] <- t[a,topo_measure]
#	}
#}
#for(a in algorithms)
#{	cat("Algorithm",a,"\n")
#	series1 <- classic_data[a,]
#	series2 <- topo_data[a,]
#	res <- t.test(x=series1, y=series2, alternative="two.sided", paired=TRUE, conf.level=0.95)
#	print(res)
#}


## display all warnings
warnings()
