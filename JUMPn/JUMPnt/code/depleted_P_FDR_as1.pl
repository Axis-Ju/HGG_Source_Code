$line=<>;
print $line;
while(<>) {
	chomp;
	@t=split /\t/,$_,-1; 
	if ($t[3]<1) {
		$t[2]=$t[11]=$t[9]=$t[10]=1;
	} 
	print join("\t",@t),"\n";
}
