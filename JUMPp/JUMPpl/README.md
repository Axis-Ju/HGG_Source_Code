# HGG_Source_Code_JUMPpl
#To localize phosphosites
We determined phosphosite reliability by the localization score (Lscore) from the JUMP software suite based on the concept of the phosphoRS algorithm.
We first derived an Lscore (0-100%) for each phosphosite in every PSM, and then aligned all phosphosites to protein sequences to produce a protein Lscore for each phosphosite. 
When one site was identified by numerous PSMs, the highest PSM Lscore was selected. 
As random assignment often occurred for ambiguous phosphosites, inflating the number of protein phosphosites, we used multiple rules to address this issue: 
(i) For ambiguous phosphosites in a PSM (e.g. the gap of the 1st and 2nd PSM Lscores < 10%), 
we searched the phosphosite information in the corresponding protein to define the site, which enabled the PSMs of low quality to borrow information from the PSMs of high quality. 
(ii) If neither the PSM Lscore nor the protein Lscore was distinguishable, we used a heuristic order to assign the phosphosite: SP-motif, S, T and Y. 
(iii) if none of the rules were applicable, we sorted the PSMs by JUMP Jscores to select the phosphosite.  
