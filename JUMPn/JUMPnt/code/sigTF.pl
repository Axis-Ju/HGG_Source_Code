$mode = $ARGV[1];
$threshold=$ARGV[2];

open(IN,"$ARGV[0]") || die "cannot open $ARGV[0]\n";
while(<IN>) {

	if (/GeneSetName/) { # header lines
		print $_;
		next;
	}

        chomp;
	$p=(split /\t/,$_)[2];
	$hit=(split /\t/,$_)[5];
	$fdr=(split /\t/,$_)[9];
	$enrichment=(split /\t/,$_)[3];

	if ($mode eq 'p') {
		if ($p<=$threshold and $enrichment>1 and $hit>1) {
			print $_,"\n";
		}
	} elsif ($mode eq 'fdr') {
		if ($fdr<=$threshold and $enrichment>1 and $hit>1) {
			print $_,"\n";
		}
	}
}
close IN;


