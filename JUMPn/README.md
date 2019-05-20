# HGG_Source_Code for JUMPn
Function of JUMPn is for network analysis of proteomics-centered multi-omics data.
JUMPn composes of 5 segments:


i. JUMPnc: Proteome co-expression clustering analysis using weighted correlation network analysis(WGCNA) published by: Peter Langfelder and Steve Horvath,WGCNA: an R package for weighted correlation network analysis.BMC Bioinformatics.2008.
R package was downlaoded: https://cran.r-project.org/web/packages/WGCNA/index.html
ii.JUMPnp: Pathway enrichment of each co-expression clusters using Fisher's exact test
iii. JUMPnk: Kinase activity reference using IKAP, published by Mischnik M et al, IKAP: A heuristic framework for inference of kinase activities from Phosphoproteomics data. Bioinformatics.2016.Downloaded at https://omictools.com/ikap-tool
iv. JUMPnt: transcriptional factor activity reference based on target gene expression
v. JUMPna: Pathway activity calculation
