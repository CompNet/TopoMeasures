TopoMeasures
=======
*Partition comparison measures for community detection*

* Copyright 2012-16 Vincent Labatut. 

TopoMeasures is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation. For source availability and license information see `licence.txt`

* Lab site: http://lia.univ-avignon.fr/
* GitHub repo: https://github.com/CompNet/TopoMeasures
* Contact: Vincent Labatut <vincent.labatut@univ-avignon.fr>

If you use this source code, please cite reference [L'15].

-----------------------------------------------------------------------


# Description
These scripts implement several measures allowing to compare two community structures, i.e. two partitions of the
node set of a given graph. They are based on popular measures defined in the field of cluster analysis, namely 
(see the source code for original bibliographic references):  

* Purity (also known under many other names in the literature, such as percent correct, accuracy, etc.)
* Rand index and its adjusted version.
* Normalized mutual information.

The scripts provided here implement (or show how to use existing implementations of) these classic measures,
and also some modified versions. Those allow to take into account some weight defined for each one of the considered elements
(in our case: nodes). The goal here is to be able to factor the network topology in these measures, which otherwise 
completely ignore this essential aspect of community detection (which is natural, since they were originally developped
to assess cluster analysis results).

The measures were originally implemented in 2012 and the corresponding paper was published in 2015 [L'15].
I did not publish the source code at this time, because I thought it was quite trivial to implement the measures I had proposed.
However, in 2016 I realized a few authors cited my work and pointed at the absence of publicly available source code, 
so I decided to clean it up and put it online. Hopefully, someone will use it! 


# Organization
The organization is very simple, all the source code is in folder `src`:

* `CommonFunctions.R`: contains some functions used to process some of the measures.
* `NormalizedMutualInformation.R`, `PurityMeasure.R` and `RandIndex.R`: process the measures listed above, as well as their modified versions. 
* `Example.R`: illustrates how to process the different measures, by applying them to
			   several data, including the example from my paper. The `Example-out-xxxxx.txt` files correspond to the results you should obtain
			   when executing this script. 


# Installation
You just need to install `R` and the required packages:

1. Install the [`R` language](https://www.r-project.org/)
2. Install the following R packages:
   * [`igraph`](http://igraph.org/r/)
3. Download this project from GitHub and unzip.
4. Launch `R`, setup the working directory with `setwd` so that it points at the root of this project. 


# Use
In order to process the measures:

1. Open the `R` console.
2. Set the project root as the working directory, using `setwd("<my directory>")`.
3. Launch the `Example.R` script, or just include one of the measure scripts (`NormalizedMutualInformation.R`, 
   `PurityMeasure.R` and `RandIndex.R`) to use in your own source code.


# Extension
The example shows how to use the measures with the three types of weights proposed in the paper [L'15].
However, the functions allow using any other scheme: you just need to process these topological weights
first, then pass them to the selected function using its `topo.measure` parameter. 


# Dependencies
* [`igraph`](http://igraph.org/r/) package: used to build and handle graphs.


# Disclaimer
This is not the exact version used in the paper [L'15], it is a reworked, cleaned and commented one.
For instance, some of the functions I had originally implemented are now available in certain libraries, 
so I decided to remove them from my scripts and use these libraries instead. 
I was careful when updating my source code, and I tested it. But of course, because of these changes, 
it is always possible that I introduced some bugs. Please, if you find one, contact me at the above
email address or post an issue on the GitHub page. 


# Reference
 * **[L'15]** V. Labatut, *Generalised measures for the evaluation of community detection methods*, International Journal of Social Network Mining (IJSNM), 2(1):44-63, 2015. [Official version](http://www.inderscienceonline.com/doi/abs/10.1504/IJSNM.2015.069776) - [arXiv version](http://arxiv.org/abs/1303.5441) - [HAL version](https://hal.archives-ouvertes.fr/hal-00802923/).
