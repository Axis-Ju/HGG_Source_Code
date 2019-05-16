
# Load the WGCNA package
#===============================================================
# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = "C:/Users/hwang1/Documents/2019_analysis/201905_HGG_Revision/whl";
setwd(workingDir); 
# Load the WGCNA package
library(WGCNA);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);





#################################    Perform WGCNA analysis  ########################################################################


# Load the WGCNA package and prepare data
library(WGCNA);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);


#### input data generated from limma.
library(dplyr)
Clus_input <- read.table("HGG_whl_input.txt",head=T,sep="\t")
head(Clus_input)
dim(Clus_input)



grep('Cortex1',colnames(Clus_input))
offst=grep('Cortex1',colnames(Clus_input))-1
offst

tb <-Clus_input[, (offst+1):(offst +10)]
dim(tb)
head(tb)

m=t(tb)
dim(m)

#  Code chunk 2
#===============================================================

gsg = goodSamplesGenes(m, verbose = 3);
gsg$allOK
dim(m)
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}


#  Code chunk 3
# enable if genes >5000
#===============================================================
enableWGCNAThreads()



#  Code chunk 
# Plot the sample tree; not nessarily needed
#===============================================================

sampleTree = hclust(dist(m), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
#par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

#ARMS and ERMS are not well seperated




#  Code chunk 5
# Choose a set of soft-thresholding powers and plot
#===============================================================
datExpr=m
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=50, by=2))
#powers = c(seq(from = 20, to=40, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType = 'signed')

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
abline(v=10,col="red")

# set the soft power
softPower = 16

#  Code chunk 5
# Evaluate adjacency and module identification
#===============================================================

# adjacency matrix
softPower
adjacency = adjacency(datExpr, power = softPower,type = "signed");

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 100;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")




#  Code chunk 6
# Calculate eigengenes and select Cutline to merge modules
#===============================================================

MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");

#select the best MEDissThres to get reasonable merged colors groups
MEDissThres = 0.2

# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;

table(mergedColors)
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

#output egigene

write.table(mergedMEs, file="hgg_whl_eigeneME_output.txt", sep="\t")

head(mergedColors)

# kME
datKME=signedKME(datExpr, mergedMEs, outputColumnName="MM.")

# attach kME to tb
nm = tb
n=ncol(nm)
nm[,(n+1):(n+ncol(datKME))]=datKME[,1:ncol(datKME)]

head(nm)


nm$group=mergedColors
head(nm)
#output clustering results
write.table(nm, file="hgg_whl_cluster_output.txt", sep="\t")

