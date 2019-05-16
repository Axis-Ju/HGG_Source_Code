

#===============================================================
# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = "C:/Users/hwang1/Documents/2016_analysis/hgg_pathway_activity_nullhypothesis testing";
setwd(workingDir);
dir()

# read log2FC between NTRK and PDGFRA
act <- read.delim("pathway_activity.txt", header = TRUE, sep = "\t")
head(act)

z <- act$LogFC

#Bootstrap to generate a list consists of 10000 acitivites each come
# from summarization of 22 randomly sampled FCs
# aP is 1.45
#Among the 22 compenents of PI3K-AKT pathway that have DE phos #change, functional annation (Ci) for protein 4,8,21 is -1, the #rest of proteins are +1
aP =list()
for (i in 1:10000) {
  a =sample(z, 22)
  a[4] <- -a[4]
  a[8] <- -a[8]
  a[21] <- -a[21]
  aPi <- sum(a)/sqrt(length(a))
  aP = c(aP, aPi)
}
length(aP)
pvalue = sum(aP > 1.45)/10000
pvalue
