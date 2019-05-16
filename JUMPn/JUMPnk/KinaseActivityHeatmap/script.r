library(limma);
#library(edgeR)
data=read.csv("IKAP_result_v1_20160422.txt", sep="\t", header=T, as.is=T);
gene=data[,1]
allDat = data;
selection = allDat;
genenames = selection[,1];
col_labels = colnames(allDat[1,]);
sampleNames = col_labels[2:length(col_labels)];
colnames(selection) = col_labels;
rownames(selection) = genenames;
A = c('PDGFRA1','PDGFRA2','PDGFRA2','PDGFRA4')
B = c('NTRK2','NTRK3')
C = c('CNTRL1','CNTRL2','CNTRL3')
mat = selection[,c(A, B, C)];
rownames(mat)=genenames
groupAOther = factor(c(rep("A", length(A)), rep("B", length(B)), rep("C", length(C))));
designA = model.matrix(~0 + groupAOther);
colnames(designA) <- c("group1", "group2", "group3");
fit <- lmFit(mat,design = designA)
contrast.matrix <- makeContrasts(group2-group1, group3-group2, group3-group1, levels=designA)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2);
options(digits=4)
top1 = topTableF(fit2,n=300000,p.value=1.0, adjust.method="BH")
passLG <- rowSums(top1[,1:3]>0.0) >= 1
top1 = top1[passLG,]
top1Pos = top1$P.Value < 0.05
top1 = top1[top1Pos,]
top1Pval = head(sort(top1$P.Value,decreasing=FALSE), n = 300000)
write.table(top1, file="Filtered.txt", sep="	");
contrast.matrix <- makeContrasts(group2-group1, group3-group2, group3-group1, levels=designA)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2);
options(digits=4)
top1 = topTableF(fit2,n=300000,p.value=1.0, adjust.method="BH")
top1Pval = head(sort(top1$P.Value,decreasing=FALSE), n = 300000)
write.table(top1, file="ALL.txt", sep="	");
library(ggplot2)

