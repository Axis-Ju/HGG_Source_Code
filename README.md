# HGG_Source_Code
The source code including two major parts: JUMPp and JUMPq.

Function of JUMPp is to process mass spectrometer generated raw data for protein identification and quantificaiton
JUMPp composes of 4 segments: 

i. Protein database search; 
ii. Putative Match filtering using a target-decoy strategy.
iii. Phosphosite site localization score calculation for phosphoproteome analysis. 
iv. protein/peptide quantification for TMT labelling

Function of JUMPn is for network analysis of proteomics-centered multi-omics data.
JUMPn composes of 5 segments:

i. Proteome co-expression clustering analysis using WGCNA: https://cran.r-project.org/web/packages/WGCNA/index.html
ii.Pathway enrichment of each co-expression clusters using Fisher's exact test
iii. Kinase activity reference using IKAP:https://omictools.com/ikap-tool
iv. transcriptional factor activity reference using target gene expression
v. Pathway activity calculation  