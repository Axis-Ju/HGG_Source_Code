java -jar DRPPM_20190514_newmachine.jar -LIMMA3Flex IKAP_result_v1_20160422.txt PDGFRA.txt NTRK_Filtered.txt CNTRL.txt Filtered.txt ALL.txt 0.05 0 false false false > script.r
R --vanilla < script.r
