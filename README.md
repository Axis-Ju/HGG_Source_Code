# HGG_Source_Code
The source code that supports the findings of this study (Deep Multiomics Profiling of Brain Tumors Identifies Signaling Networks Downstream of Cancer Driver Genes) includes two major parts: JUMPp and JUMPq.

Function of JUMPp is to pre-process mass spectrometer generated raw data for protein identification and quantificaiton
JUMPp composes of 4 segments: 

i. JUMPps: Protein database search; 
ii. JUMPpf: Putative Match filtering using a target-decoy strategy.
iii. JUMPpl: Phosphosite site localization score calculation for phosphoproteome analysis. 
iv. JUMPpg: protein/peptide quantification for TMT labelling


Function of JUMPn is for network analyses of post-processed proteomics-centered multi-omics data.
JUMPn composes of 5 segments:

i. JUMPnc: Proteome co-expression clustering analysis using weighted correlation network analysis(WGCNA) published by: Peter Langfelder and Steve Horvath,WGCNA: an R package for weighted correlation network analysis.BMC Bioinformatics.2008.
R package was downlaoded: https://cran.r-project.org/web/packages/WGCNA/index.html
ii.JUMPnp: Pathway enrichment of each co-expression clusters using Fisher's exact test
iii. JUMPnk: Kinase activity reference using IKAP, published by Mischnik M et al, IKAP: A heuristic framework for inference of kinase activities from Phosphoproteomics data. Bioinformatics.2016.Downloaded at https://omictools.com/ikap-tool
iv. JUMPnt: transcriptional factor activity reference based on target gene expression
v. JUMPna: Pathway activity calculation

JUMP softwares that produced by Peng lab were developed by: Xusheng Wang, Ji-Hoon Cho, Yuxin Li, Timothy Shaw, Hong Wang and Junmin Peng  