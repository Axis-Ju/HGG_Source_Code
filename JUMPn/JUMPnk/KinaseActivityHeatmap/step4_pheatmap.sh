java -jar DRPPM_20190514_newmachine.jar -ExtractDEGenes Filtered.txt > anova_geneList.txt

java -jar DRPPM_20190514_newmachine.jar -FilterKinaseBasedOnFrequency ../IKAP_Input_20160420/jpeng_phosphosite_matrix_normalized_v1_20160420.txt anova_geneList.txt > anova_geneList_3sub.txt

java -jar DRPPM_20190514_newmachine.jar -plotPHeatMap IKAP_result_v1_20160422.txt sampleNames.txt anova_geneList_3sub.txt IKAP_anova_v1_20160705.png Kinase_Activity true false false 500 800 12 12 > pheatmap.r
R --vanilla < pheatmap.r

# note we edited the pheatmap.r to add EEF2k and ULK1 the final version is final_pheatmap.r

