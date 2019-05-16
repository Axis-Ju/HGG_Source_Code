allDat = read.table("IKAP_result_v1_20160422.txt", header=TRUE, row.names=1 );
selection = allDat;
A = c('CNTRL1','CNTRL2','CNTRL3','PDGFRA1','PDGFRA2','PDGFRA3','PDGFRA4','NTRK1','NTRK2','NTRK3')
col_labels = A;
sampleLocation = col_labels;
sampleNames = col_labels;
labels = sampleNames;
colnames(selection) = col_labels;
rows = length(selection[,1]);
zselection = apply(selection, 1, scale);
zselection = apply(zselection, 1, rev)
colnames(zselection) = names(selection)
selection = as.matrix(zselection);
geneList = "MAPK7";
geneList = c(geneList, "PRKDC");
geneList = c(geneList, "PRKCZ");
geneList = c(geneList, "PRKCI");
geneList = c(geneList, "PRKAA2");
geneList = c(geneList, "PRKAA1");
geneList = c(geneList, "SGK1");
geneList = c(geneList, "MAP2K7");
geneList = c(geneList, "CDK5");
geneList = c(geneList, "CHEK1");
geneList = c(geneList, "CAMK2A");
geneList = c(geneList, "PRKACA");
geneList = c(geneList, "MAPK3");
geneList = c(geneList, "MAP2K1");
geneList = c(geneList, "NUAK1");
geneList = c(geneList, "PRKCE");
geneList = c(geneList, "DYRK2");
geneList = c(geneList, "PRKCG");
geneList = c(geneList, "MAPK11");
geneList = c(geneList, "PAK3");
geneList = c(geneList, "SRC");
geneList = c(geneList, "UHMK1");
geneList = c(geneList, "AURKB");
geneList = c(geneList, "MAPK13");
geneList = c(geneList, "PRKCA");
geneList = c(geneList, "AKT1");
geneList = c(geneList, "PRKG2");
geneList = c(geneList, "MAP2K4");
geneList = c(geneList, "CDK2");
geneList = c(geneList, "FYN");
geneList = c(geneList, "ATR");
geneList = c(geneList, "CHEK2");
geneList = c(geneList, "CAMK2D");
geneList = c(geneList, "GSK3A");
geneList = c(geneList, "ATM");
geneList = c(geneList, "MAPKAPK2");
geneList = c(geneList, "PAK1");
geneList = c(geneList, "RPS6KB1");
geneList = c(geneList, "CSNK1D");
geneList = c(geneList, "EEF2K");
geneList = c(geneList, "ULK1");
dataset = selection[geneList,]
library(pheatmap)
minimum = -3;
maximum = 3;
if (abs(min(dataset)) > abs(max(dataset))) {
dataset[dataset < -abs(max(dataset))] = -abs(max(dataset))
} else {
dataset[dataset > abs(min(dataset))] = abs(min(dataset))
}
bk = c(seq(minimum,minimum/2, length=100), seq(minimum/2,maximum/2,length=100),seq(maximum/2,maximum,length=100))
hmcols<- colorRampPalette(c("#34c5fd","black","red"))(length(bk)-1)
png(file = "IKAP_anova_v1_20160705.png", width=500,height=800)
pheatmap(dataset, cluster_col = F, cluster_row = T, fontsize_row = 12, fontsize_col = 12, show_rownames = T, color=hmcols)
dev.off();

