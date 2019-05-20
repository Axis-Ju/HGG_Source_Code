#HGG_Source_Code_TF activity reference.
The Purpose of the code is to reference TF activity based on target gene expression.
Required to install 'drppm' package before usage:  https://github.com/gatechatl/DRPPM

Usage:
	perl jump_nenrich.pl -R_code_path filterFDRmatrix.R -foreground <foreground-genes> -background <background-genes> -gene-set <pathway-DB>

Example: 
	cp example_data/gene_list_clustered.txt .
	perl code/jump_nenrich.pl -R_code_path code/filterFDRmatrix.R -foreground gene_list_clustered.txt -background example_data/hs_mm_homo_r66.txt -gene-set TF_DB/TF_target_list.txt

