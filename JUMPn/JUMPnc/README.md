# HGG_Source_Code_JUMPnc 
The DE proteome and phosphoproteome were analyzed by WGCNA co-expression clustering analysis. 
Briefly, Pearson correlation matrix was calculated using all samples, allowing only positive correlation. 
Hybrid dynamic tree-cutting method with a minimum height for merging modules at 0.15 was applied to define co-expression clusters. 
The first principal component (i.e. eigengene) was calculated as a consensus trend for each co-expression cluster. DE Proteins were assigned to each cluster based on Pearson r values. 
WGCNA was published by: Langfelder P, Horvath S. WGCNA: an R package for weighted correlation network analysis. BMC bioinformatics 9, 559 (2008)
The source code for functions necessary to perform Weighted Correlation Network Analysis was published:https://cran.r-project.org/web/packages/WGCNA/index.html