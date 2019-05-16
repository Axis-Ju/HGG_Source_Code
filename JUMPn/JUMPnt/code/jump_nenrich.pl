#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;

use FindBin qw($Bin);
use lib "$Bin";

my @params = (
"-R_code_path=s",
"-foreground=s",
"-background=s",
"-gene-set=s", # e.g., GO database; could be multiple files and comma separated
"-FDR-cutoff=s", # defualt = 0.2, i.e., 20% BH FDR
"-P-cutoff=s", # defualt = 0.05
"-rank-cutoff=s", # top how many GO terms are shown in the heatmap
"-filter-method=s", # 'rank' (actually FDR + rank) or 'FDR' (only FDR; good for TF target analysis)
);

my (%FLAGS,%parahash);
GetOptions(\%FLAGS, @params);

$parahash{'code_path'}=$Bin;

if ($FLAGS{"R_code_path"}){
        $parahash{"R_code_path"}=$FLAGS{"R_code_path"};
}else {
        die "specify -R_code_path [absolute path to filterFDRmatrix.R]\n";;
}

if ($FLAGS{"filter-method"}){
        $parahash{"filter-method"}=$FLAGS{"filter-method"};
}else {
        $parahash{"filter-method"}='rank';
}

if ($FLAGS{"rank-cutoff"}){
        $parahash{"rank-cutoff"}=$FLAGS{"rank-cutoff"};
}else {
        $parahash{"rank-cutoff"}=10;
}

if ($FLAGS{"FDR-cutoff"}){
        $parahash{"FDR-cutoff"}=$FLAGS{"FDR-cutoff"};
}else {
        $parahash{"FDR-cutoff"}=0.20;
}

$parahash{"P-cutoff"}=(defined($FLAGS{"P-cutoff"}))?$FLAGS{"P-cutoff"}:0.05;

if ($FLAGS{"foreground"}){
        $parahash{"foreground"}=$FLAGS{"foreground"};
}else {
        die "specify -foreground [foreground gene list file]\n";
}

if ($FLAGS{"background"}){
        $parahash{"background"}=$FLAGS{"background"};
}else {
        die "specify -background [background gene list file]\n";
}

if ($FLAGS{"gene-set"}){
        $parahash{"gene-set"}=$FLAGS{"gene-set"};
}else {
        die "specify -gene-set [gene-set gene list file]\n";
}

# multiple annotation database files:
my @annotationDB=split /\,/,$parahash{"gene-set"};

for (my $i=0; $i<=$#annotationDB; $i++) {
	print "Using the following database: $annotationDB[$i]\n";
	# Fisher.s exact test
	my $cmd="drppm -OverRepresentationAnalysisWithoutFilter $parahash{foreground} $annotationDB[$i] $parahash{background} false $parahash{foreground}_vs_DB$i\_enrichAnalysis.txt";
	print "$cmd\n";
	system(qq(drppm -OverRepresentationAnalysisWithoutFilter $parahash{"foreground"} $annotationDB[$i] $parahash{"background"} false $parahash{"foreground"}_vs_DB$i\_enrichAnalysis.txt >/dev/null 2>&1));
	# FDR
	system(qq(drppm -OverRepresentationAnalysisFDR $parahash{"foreground"}_vs_DB$i\_enrichAnalysis.txt $annotationDB[$i] $parahash{"foreground"}_vs_DB$i\_enrichAnalysis_FDR.txt >/dev/null 2>&1));
}
system(qq(cat $parahash{"foreground"}_vs_DB*enrichAnalysis_FDR.txt > $parahash{"foreground"}_enrichAnalysis_FDR.txt));

# Filter list by BH FDR
#print "test\n";
system(qq(perl $parahash{'code_path'}/sigTF.pl $parahash{"foreground"}_enrichAnalysis_FDR.txt fdr $parahash{"FDR-cutoff"} > $parahash{"foreground"}_enrichAnalysis_FDRfiltered.txt));
system(qq(perl $parahash{'code_path'}/sigTF.pl $parahash{"foreground"}_enrichAnalysis_FDR.txt p $parahash{"P-cutoff"} > $parahash{"foreground"}_enrichAnalysis_Pfiltered.txt));

# prepare publication table:
# 1) separate different DBs: if header row, print "\n",$_;
# 2) re-organize columns
publicationTable("$parahash{foreground}_enrichAnalysis_FDRfiltered.txt","$parahash{foreground}_Publication_table_FDRfiltered.txt");
publicationTable("$parahash{foreground}_enrichAnalysis_Pfiltered.txt","$parahash{foreground}_Publication_table_Pfiltered.txt");

# depleted -> p/FDR = 1
system(qq(perl $parahash{'code_path'}/depleted_P_FDR_as1.pl $parahash{"foreground"}_enrichAnalysis_FDR.txt > $parahash{"foreground"}_enrichAnalysis_FDR_smDpl.txt));
publicationTable("$parahash{foreground}_enrichAnalysis_FDR_smDpl.txt","$parahash{foreground}_Publication_table.txt");
# filter by p value
#system(qq(cat $parahash{"foreground"}_enrichAnalysis_FDR_smDpl.txt | perl -lane \'chomp;\@t=split /\\t/;if ($t[2] < $parahash{"P-cutoff"})));
=head
# prepare FDR matrix
system(qq(mkdir FDRmatrix; cp $parahash{"foreground"}_enrichAnalysis_FDR.txt FDRmatrix));
# split files
system(qq(cd FDRmatrix; perl $parahash{'code_path'}/splitFile.pl $parahash{"foreground"}_enrichAnalysis_FDR.txt));
# edil column names
system(qq(cd FDRmatrix; ls *enrichAnalysis.txt > input.txt));
# edit column names
system(qq(cd FDRmatrix; perl $parahash{'code_path'}/editColnames.pl input.txt));
# combine tables
system(qq(cd FDRmatrix; perl $parahash{'code_path'}/combineTables.pl *updated.txt));

# Filter FDR matrix, print heatmap
# rm '#' in header
system(qq(cd FDRmatrix; cat combineTable_for_heatmap.txt | perl -lane 'chomp;s/\#//g;print;' > combineTable_for_heatmap2.txt));
# R:
system(qq(cd FDRmatrix; R CMD BATCH --no-save --args -method=$parahash{"filter-method"} -combineTable_for_heatmap=combineTable_for_heatmap2.txt -FDRcutoff=$parahash{"FDR-cutoff"} -rankcutoff=$parahash{"rank-cutoff"} $parahash{'R_code_path'} R.log));
=cut
#----------------------------------------------------------------------
sub publicationTable {
# prepare publication table:
# 1) separate different DBs: if header row, print "\n",$_;
# 2) re-organize columns
	my ($input,$output)=@_;

	open(IN,$input) || die "cannot open $input\n";
	open(OUT,"\>",$output);

	while(<IN>) {
		chomp;
		my @t=split /\t/,$_;
		if ($t[1] eq 'GeneSetName') {
			print OUT "\n";
		}
		print OUT "$t[0]\t$t[1]\t$t[6]\t$t[7]\t$t[5]\t$t[4]\t$t[3]\t$t[2]\t$t[9]\n";
	}
	close IN;
	close OUT;
}
