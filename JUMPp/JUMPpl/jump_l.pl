#!/usr/bin/perl

use Getopt::Long;

use Cwd;
use Cwd 'abs_path';
use File::Basename;
use FindBin qw($Bin);
use lib "$Bin";
use Clone qw(clone);
use Spiders::ProcessingMzXML;
use Spiders::Params;
use List::Util qw(shuffle);
use Parallel::ForkManager;

my ($help,$parameter,$raw_file);
GetOptions('-help|h'=>\$help,
	'-p=s'=>\$parameter,
);

my $progname = $0;
# remove path from our name
$progname =~ s@(.*)/@@i;

usage() if ($help || !defined($parameter));
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time); 
my $Hydrogen_mass = 1.007825032;

##################################################
## Read a parameter file and initialization	##
##################################################
my $p = Spiders::Params->new('-path'=>$parameter);
my $params = $p->parse_param();
$params->{'Mn_Mn1'} = 0.5;
$params->{'M_M1'} = 0.3;
$params->{'low_mass'} = 57;
$params->{'high_mass'} = 187;
my $IDmod = $params->{'IDmod'};
my $Outfile = $params->{'Output'};
my $dir = (split(/\//, $IDmod))[-2];
$dir =~ s/^sum/loc/;
$dir = getcwd()."/".$dir;
if (!-e $dir) {
	system(qq(mkdir $dir >/dev/null 2>&1));
}
system(qq(cp $parameter $dir));
my $JUMPLRES = $dir."/.ID.lscore";
open(JUMPLRES,">", $JUMPLRES);
my $logFile = $dir."/jump_l.log";
open(LOGFILE,">", $logFile);
print "\n\n  Start: ";
printf "%4d-%02d-%02d %02d:%02d:%02d\n", $year + 1900, $mon + 1, $mday, $hour, $min, $sec;
print LOGFILE "\n\n  Start: ";
printf LOGFILE "%4d-%02d-%02d %02d:%02d:%02d\n", $year + 1900, $mon + 1, $mday, $hour, $min, $sec;
print "  Initializing jump -l program\n\n";
print LOGFILE "  Initializing jump -l program\n\n";

##########################
## Read IDmod.txt file	##
##########################
my ($frac_scan,$databaseHeader) = read_IDmod($IDmod);
print "\n";
print LOGFILE "\n";

##########################################################################################
## For each fraction, extract scan information from a .mzXML file and create .dta files	##
##########################################################################################
my $MAX_PROCESSES = 32;
print "  Extracting m/z and generating dta files for each fraction \n";
print LOGFILE "  Extracting m/z and generating dta files for each fraction \n";
if ($params -> {'cluster'} eq '1') {
	my $nTotalJobs = 0;
	my %jobIDs;
	foreach my $fraction (keys %$frac_scan) {
		my $basename = basename($fraction);
		my $new_path = $dir."/$basename";
		system(qq(mkdir $new_path >/dev/null 2>&1));
		system(qq(cp $parameter "$new_path/" >/dev/null 2>&1));	
		my $jobName = "Extraction_".$nTotalJobs;
		open (JOB, ">", "$new_path/$jobName.sh") or die "Cannot creat a job file\n";
		print JOB "#!/bin/bash\n";
		print JOB "#\$ -N $jobName\n";
		print JOB "#\$ -e $new_path/$jobName.e\n";
		print JOB "#\$ -o $new_path/$jobName.o\n";
		print JOB "perl $Bin/Extraction_runshell.pl -fraction $fraction -outdir $dir -IDmod $IDmod -parameter $parameter\n\n";
		close (JOB);
		my $command = qq(qsub -cwd -pe mpi 8 -l mem_free=16G,h_vmem=16G "$new_path/$jobName.sh");
		my $job = lc(qx[$command]);
		chomp ($job);
		if ($job =~ /job (\d+)/) {
			$jobIDs{$1} = 1;
		}
		$nTotalJobs++;
		print "\r  $nTotalJobs jobs are submitted";
	}
	print "\n  You submitted $nTotalJobs job(s) for data extraction \n";
	print LOGFILE "  You submitted $nTotalJobs job(s) for data extraction \n";
	CheckJobStat($nTotalJobs, \%jobIDs);
} elsif ($params -> {'cluster'} eq '0') {
	my $nTotalJobs = 0;
	my $pm = new Parallel::ForkManager($MAX_PROCESSES);
	foreach my $fraction (keys %$frac_scan) {	
		$nTotalJobs++;
		print "\r  $nTotalJobs job(s) submitted";
		$pm -> start and next;
		my $basename = basename($fraction);
		my $new_path = $dir."/$basename";
		system(qq(mkdir $new_path >/dev/null 2>&1));
		system(qq(cp $parameter "$new_path/" >/dev/null 2>&1));
		my $jobName = "Extraction_" . $nTotalJobs;
		open (JOB, ">", "$new_path/$jobName.sh") or die "Cannot creat a job file\n";
		print JOB "#!/bin/bash\n";
		print JOB "#\$ -N $jobName\n";
		print JOB "#\$ -e $new_path/$jobName.e\n";
		print JOB "#\$ -o $new_path/$jobName.o\n";
		print JOB "perl $Bin/Extraction_runshell.pl -fraction $fraction -outdir $dir -IDmod $IDmod -parameter $parameter\n\n";
		close (JOB);
		system(qq(sh "$new_path/$jobName.sh" > /dev/null 2>&1));
		$pm -> finish;
	}
	$pm -> wait_all_children;
	print "\n  $nTotalJobs job(s) finished\n";
	print LOGFILE "\n  $nTotalJobs job(s) finished\n";
}

##########################################################################################
## For each peptide, calculate localization scores of all possible modified peptides	##
##########################################################################################
print "\n\n  Calculating localization scores\n";
print LOGFILE "\n\n  Calculating localization scores";
if ($params -> {'cluster'} eq '1') {
	my $nEntriesPerJob = 100;
	my $nTotalJobs = 0;
	my $maxJobs = 100;
	my %jobIDs;
	foreach my $frac (keys %$frac_scan) {
		my $basename = basename($frac);
		my $new_path = $dir."/$basename";
		my @outfiles = keys %{$frac_scan->{$frac}};
		my $nEntries = scalar(@outfiles);
		my $nJobs = int($nEntries / $nEntriesPerJob) + 1;
		if ($nJobs > $maxJobs) {
			$nEntriesPerJob = int($nEntries / $maxJobs) + 1;
			$nJobs = int($nEntries / $nEntriesPerJob) + 1;
		}      
		for (my $i = 0; $i < $nJobs; $i++) {
			my $jobName = "Job_PTM_".$nTotalJobs;
			open (JOB, ">", "$new_path/$jobName.sh") or die "Cannot creat a job file\n";
			print JOB "#!/bin/bash\n";
			print JOB "#\$ -N $jobName\n";
			print JOB "#\$ -e $new_path/$jobName.e\n";
			print JOB "#\$ -o $new_path/$jobName.o\n";
			for (my $j = 0; $j < $nEntriesPerJob; $j++) {
				my $k = $nEntriesPerJob * $i + $j;
				last if ($k >= $nEntries);
				my $queryOutfile =  $outfiles[$k];
				my $queryPeptide = "\"" . $frac_scan->{$frac}->{$queryOutfile}->{'peptide'} . "\"";
				print JOB "perl /data1/pipeline/release/version12.1.0/JUMPl/JUMPl_runshell.pl -fraction $frac -outdir $new_path -scan $queryOutfile -peptide $queryPeptide -parameter $parameter\n\n";
			}
			close (JOB);
			my $command = qq(qsub -cwd "$new_path/$jobName.sh");
		        my $job = lc(qx[$command]);
	        	chomp ($job);
		        if ($job =~ /job (\d+)/) {
	        	        $jobIDs{$1} = 1;
		        } 
			$nTotalJobs++;
			print "\r  $nTotalJobs jobs are submitted";
		}
	}
	print "\n  You submitted $nTotalJobs job(s) for local scoring \n";
	print LOGFILE "\n  You submitted $nTotalJobs job(s) for local scoring \n";
	CheckJobStat($nTotalJobs, \%jobIDs);
} elsif ($params -> {'cluster'} eq '0') {
	my $nEntriesPerJob = 100;
	my $nTotalJobs = 0;
	my $maxJobs = 100;
	$pm = new Parallel::ForkManager($MAX_PROCESSES);
	foreach my $frac (keys %$frac_scan) {
		my $basename = basename($frac);
		my $new_path = $dir."/$basename";
		my @outfiles = keys %{$frac_scan->{$frac}};
		my $nEntries = scalar(@outfiles);
		my $nJobs = int($nEntries / $nEntriesPerJob) + 1;
		if ($nJobs > $maxJobs) {
			$nEntriesPerJob = int($nEntries / $maxJobs) + 1;
			$nJobs = int($nEntries / $nEntriesPerJob) + 1;
		}
		for (my $i = 0; $i < $nJobs; $i++) {
			$nTotalJobs++;
			print "\r  $nTotalJobs job(s) submitted (some of them may be finished)";
			my $jobName = "Job_PTM_".$nTotalJobs;
			$pm -> start and next;
			open (JOB, ">", "$new_path/$jobName.sh") or die "Cannot creat a job file\n";
			print JOB "#!/bin/bash\n";
			print JOB "#\$ -N $jobName\n";
			print JOB "#\$ -e $new_path/$jobName.e\n";
			print JOB "#\$ -o $new_path/$jobName.o\n";
			for (my $j = 0; $j < $nEntriesPerJob; $j++) {
				my $k = $nEntriesPerJob * $i + $j;
				last if ($k >= $nEntries);
				my $queryOutfile =  $outfiles[$k];
				my $queryPeptide = "\"" . $frac_scan->{$frac}->{$queryOutfile}->{'peptide'} . "\"";
				print JOB "perl /data1/pipeline/release/version12.1.0/JUMPl/JUMPl_runshell.pl -fraction $frac -outdir $new_path -scan $queryOutfile -peptide $queryPeptide -parameter $parameter\n\n";
			}
			close (JOB);
			system(qq(sh "$new_path/$jobName.sh" > /dev/null 2>&1));
			$pm -> finish;
		}		
	}
	$pm -> wait_all_children;
	print "\n  $nTotalJobs job(s) finished\n";
	print LOGFILE "\n  $nTotalJobs job(s) finished\n";
}

sub CheckJobStat {
	my ($nJobs, $jobIDs) = @_;
	my $nFinishedJobs = 0;
	my $jobInfo = 1;
	my ($username) = getpwuid($<);
	while($jobInfo) {
		my $command =  "qstat -u $username";
		my $jobStatus = qx[$command];
		my @jobStatusArray = split(/\n/,$jobStatus);
		if (@jobStatusArray) {	# i.e. There are some jobs in the queue (may include other jobs like searching)
			my @jobIDsInQueue;
			foreach (@jobStatusArray) {
				my ($jobID) = ($_ =~ /([0-9]+)\s/);
				if (defined $$jobIDs{$jobID}) {
					push (@jobIDsInQueue, $jobID);
				}
			}
			if (@jobIDsInQueue) { # i.e. There are jobs of interest in the queue
				$nFinishedJobs = $nJobs - scalar(@jobIDsInQueue);
				print "\r  $nFinishedJobs jobs are finished";
				if ($nFinishedJobs == $nJobs) {
					$jobInfo = 0;
				} else {
					$jobInfo = 1;
				}
			} else {	# i.e. There are jobs in the queue, but all jobs of interest are finished
				print "\r  $nJobs jobs are finished";
				$jobInfo = 0;
			}
		} else {	# i.e. There's no job in the queue
			print "\r  $nJobs jobs are finished";
			$jobInfo = 0;
		}
		sleep(5);
	}
}

##########################################################
## Create .ID.lscore file (a hidden intermediate file)	##
##########################################################
print "\n";
print JUMPLRES "Peptide;Protein;Outfile;measuredMH;calcMH;ppm;XCorr;dCn;Ions;red;group;subgroup;unique;tryptic;pos;precursor_peak_intensity_percentage;JUMPl_peptide;JUMPl_site;JUMPl_score\n";
foreach my $fraction (keys %$frac_scan) {
	my $basename = basename($fraction);
	my $new_path = $dir."/$basename";
	foreach my $scan (keys %{$frac_scan->{$fraction}}) {
		print "\r  Summarizing scans: $scan";
		open(INFILE,"$new_path/$scan.out");
		my $data = <INFILE>;
		chomp $data;
		my @data_array = split(/\t/,$data);

		foreach my $line (@{$frac_scan->{$fraction}->{$scan}->{'ALL'}}) {
			my ($AAstartPosition) = (split(/\;/, $line))[14] =~ /^AA(\d+)to/;
			my $Jumpl_result = join(';',@data_array[3..$#data_array]);
			print JUMPLRES $line,";",$Jumpl_result,"\n";
		}
	}
}
print "\n\n";
print LOGFILE "\n\n";
close(JUMPLRES);

##########################################################
## Retrieve the information from .ID.lscore file and	##
## create hashes containing modification site scores	##
##########################################################
my $inputFile = $JUMPLRES;
my $modSymbols = "#%*";
my $tol = $params->{'pertide_score_tolerance'};
open (IDmod, "<", $inputFile) or die "Cannot open $inputFile\n";
my $header = <IDmod>;
my %pepScores;
my %protScores;

print "  Retrieving information from the summarized file (may take a while)\n";
print LOGFILE "  Retrieving information from the summarized file (may take a while)\n";
my $nLines = 0;
while (<IDmod>) {
	$nLines++;
	print "\r  Gathering information from $nLines PSMs";
	chomp;
	my $line = $_;
	my @elems = split(/;/, $line);
	my $peptide = $elems[0];
	my $protein = $elems[1];
	my $outfile = $elems[2];
	
	## NOTE: Assumed that modification site information starts from the 17th column of a table
	## Check function for the existence of modification site information
	if (!defined $elems[16]) {
		print "  No modification site information\n";
		print "  Check the following entry\n";
		print "  $line\n";
		print LOGFILE "  No modification site information\n";
		print LOGFILE "  Check the following entry\n";
		print LOGFILE "  $line\n";
		exit;
	}

	## Generate hashes containing site scores at the peptide and protein levels
	$pepScores{$peptide}{$outfile}{'info'} = $line;
	for (my $i = 16; $i < scalar(@elems); $i = $i + 2) {
		## Find the site location in the protein
		my ($AAStartPos) = ($elems[14] =~ /AA(\d+)to/);
		my ($modAA, $pepPos) = ($elems[$i] =~ /([A-Z])(\d+)/);
		my $protSite = $modAA.($AAStartPos + $pepPos - 1);
		## Site score(s) at the peptide level
		$pepScores{$peptide}{$outfile}{'sites'}{$elems[$i]} = $elems[$i + 1];
		$pepScores{$peptide}{$outfile}{'proteins'}{$protein}{$elems[$i]} = $protSite; 
		## Site score(s) at the protein level
		if (!defined $protScores{$protein}{$protSite}) {
			$protScores{$protein}{$protSite} = $elems[$i + 1];
		} else {
			if ($elems[$i + 1] > $protScores{$protein}{$protSite}) {
				# Update the site score at the protein level
				# Take the highest site score value from the corresponding peptide scores
				$protScores{$protein}{$protSite} = $elems[$i + 1];
			}
		}
	}
}
close (IDmod);
print "\n";
print LOGFILE "  Gathered information from $nLines PSMs\n";

##########################################################
## Update the modification site(s) based on site scores	##
##########################################################
my $nRandom = 0;
my $nTotPeptides = scalar(keys %pepScores);
my $nPeptides = 0;
foreach my $peptide (keys %pepScores) {
	$nPeptides++;
	print "\r  Updating the modification-site scores in $nPeptides peptides (out of $nTotPeptides peptides)";

	## 1. Count the number and position(s) of modification site(s) in the original peptide
	my $origPeptide = $peptide;
	$origPeptide = (split(/\./, $origPeptide))[1];
	$origPeptide =~ s/\@//g;	# Remove methionine modification symbol (@)
	my $nOrigSites = 0;
	my @origSites;
	while ($origPeptide =~ /[$modSymbols]/g) {
		## n : the number of modification sites in the original peptide (=$nSites)
		$nOrigSites++;
		my $pos = pos ($origPeptide) - $nOrigSites;
		my $AA = substr($origPeptide, $pos + $nOrigSites - 2, 1);
		push (@origSites, $AA.$pos);
	}

	## 2. Go into each PSM and update the modification site(s)
	while (my ($outfile, $outfileHash) = each %{$pepScores{$peptide}}) {
		my $nSites = $nOrigSites;
		my $nSTYs = scalar(keys %{$$outfileHash{'sites'}});
		## When there's no possible modification sites other than the identified one(s),
		## we can use the original site(s)/score(s) as a localization result
		## e.g. peptide = ABCDES#FGHI
		if ($nSTYs == $nSites) {
			$$outfileHash{'selectedSites'} = \%{$$outfileHash{'sites'}};
			$$outfileHash{'selectedLevel'} = "peptide";
			next;
		}
		## When there are possible modification sites other than the identified one(s)
		## e.g. peptide = ABCDS#EFGSHITJKL
		my @pepSiteScores;
		foreach my $site (keys %{$$outfileHash{'sites'}}) {
			push (@pepSiteScores, $$outfileHash{'sites'}{$site});
		}
		@pepSiteScores = sort {$b <=> $a} @pepSiteScores;

		##########################################################################
		## 1. Choose sites as many as possible based on peptide-level scores	##
		##########################################################################
		## Select possible site-scores from @pepSiteScores
		for (my $i = $nSites - 1; $i < scalar(@pepSiteScores) - 1; $i++) {
			if ($pepSiteScores[$i] > $pepSiteScores[$i + 1] + $tol) {
				@pepSiteScores = @pepSiteScores[0..$i];
				last;
			}
		}
		## Clearly, all sites can be clearly selected based on peptide-level scores
		my $nSelected = 0;
		if (scalar(@pepSiteScores) == $nSites) {
			for (my $i = 0; $i < $nSites; $i++) {
				foreach my $site (keys %{$$outfileHash{'sites'}}) {
					if ($$outfileHash{'sites'}{$site} == $pepSiteScores[$i]) {
						$$outfileHash{'selectedSites'}{$site} = $pepSiteScores[$i];
						$$outfileHash{'selectedLevel'} = "peptide";
						$nSelected++;
					}
				}
			}
			next;
		## Some sites are ambiguous in terms of peptide-level scores
		## Select as many sites as possible at the peptide level
		} elsif (scalar(@pepSiteScores) > $nSites) {
			for (my $i = 0; $i < $nSites; $i++) {
				if ($pepSiteScores[$i] > $pepSiteScores[$nSites] + $tol) {
					foreach my $site (keys %{$$outfileHash{'sites'}}) {
						if ($$outfileHash{'sites'}{$site} == $pepSiteScores[$i]) {
							next if (defined $$outfileHash{'selectedSites'}{$site});
							$$outfileHash{'selectedSites'}{$site} = $pepSiteScores[$i];
							$nSelected++;
						}
					}
				}
			}
		} else {
			print "  Fewer peptide-level site scores than the modification sites\n";
			print "  Check $peptide in $outfile\n";
			print LOGFILE "  Fewer peptide-level site scores than the modification sites\n";
			print LOGFILE "  Check $peptide in $outfile\n";
			exit;
		}
		if ($nSelected == $nSites) {
			## Exit if all sites are determined based on peptide-level scores
			next;
		} elsif ($nSelected > 0 && $nSelected < $nSites) {
			## If some sites are determined based on peptide-level scores, those sites should not be considered hereafter
			$nSites = $nSites - $nSelected;
			for (my $i = 0; $i < $nSelected; $i++) {
				shift @pepSiteScores;
			}
			## Go to the next step; consider protein-level site scores and choose modification sites
		} elsif ($nSelected > $nSites) {
			print "  More peptide-level site scores are selected than the modification sites\n";
			print "  Check $peptide in $outfile\n";
			print LOGFILE "  More peptide-level site scores are selected than the modification sites\n";
			print LOGFILE "  Check $peptide in $outfile\n";
			exit;
		}

		##################################################################################
		## 2. Since some sites have the same peptide-level scores, 			##
		##    introduce protein-level scores in order to choose modification sites	##
		##################################################################################
		## Find the peptide modification sites to be considered
		## In @pepSiteScores, very low scored sites are already eliminated
		## 1 modification site: sites with scores >= (highest score - tolerance) will be considered hereafter
		## 2 modification sites: sites with scores >= (second highest score - tolerance) will be considered hereafter
		## 3 modification sites: sites with scores >= (third highest score - tolerance) will be considered hereafter
		## etc.
		my @pepPossibleSites;
		foreach my $site (keys %{$$outfileHash{'sites'}}) {
			if ($$outfileHash{'sites'}{$site} >= $pepSiteScores[$nSites - 1] - $tol) {
				next if (defined $$outfileHash{'selectedSites'}{$site});
				push (@pepPossibleSites, $site);
			}
		}
		## Protein-level scores of those modification sites
		my %protSiteScores;
		foreach my $protein (keys %{$$outfileHash{'proteins'}}) {
			foreach my $pepPossibleSite (@pepPossibleSites) {
				my $protSite = $$outfileHash{'proteins'}{$protein}{$pepPossibleSite};
				if (!defined $protSiteScores{$pepPossibleSite}) {
					$protSiteScores{$pepPossibleSite} = $protScores{$protein}{$protSite};
				} else {
					if ($protSiteScores{$pepPossibleSite} < $protScores{$protein}{$protSite}) {
						$protSiteScores{$pepPossibleSite} = $protScores{$protein}{$protSite};
					}
				}
			}
		}
		my @protSiteScoresArray;
		foreach my $site (keys %protSiteScores) {
			push (@protSiteScoresArray, $protSiteScores{$site});
		}
		@protSiteScoresArray = sort {$b <=> $a} @protSiteScoresArray;
		
		## Choose sites whose protein-level scores are much greater than others
		$nSelected = 0;
		for (my $i = 0; $i < $nSites; $i++) {
			if ($protSiteScoresArray[$i] > $protSiteScoresArray[$i + 1]) {
				foreach my $site (keys %protSiteScores) {
					if ($protSiteScores{$site} == $protSiteScoresArray[$i]) {
						$$outfileHash{'selectedSites'}{$site} = $$outfileHash{'sites'}{$site};
						$nSelected++;
						## If some sites are determined based on protein-level scores, those sites should not be considered hereafter
						delete $protSiteScores{$site};
					}
				}
			}
		}
		if ($nSelected == $nSites) {
			## If all sites are determined based on protein-level scores
			## 1. Assign alternative sites
			## 2. Go to the next outfile or peptide
			foreach my $site (keys %protSiteScores) {
				$$outfileHash{'alternativeSites'}{$site} = $$outfileHash{'sites'}{$site};
			}
			next;
		} elsif ($nSelected > 0 && $nSelected < $nSites) {
			## If some sites are determined based on protein-level scores, those sites should not be considered hereafter
			$nSites = $nSites - $nSelected;
			for (my $i = 0; $i < $nSelected; $i++) {
				shift @protSiteScoresArray;
			}
			## Go to the next step; look for "SP" sites and choose them if exist
		} elsif ($nSelected > $nSites) {
                        print "  More peptide- and protein-level site scores are selected than the modification sites\n";
                        print "  Check $peptide in $outfile\n";
                        print LOGFILE "  More peptide- and protein-level site scores are selected than the modification sites\n";
                        print LOGFILE "  Check $peptide in $outfile\n";
			exit;
		}

		##################################################################################
		## 3. Some sites still have the same scores even at the protein-level			##
		##    We have to apply our own rules to select reasonable modification sites	##
		##################################################################################
		## Rule 1. Choose "SP" sites
		$nSelected = 0;
		my @SPsites;
		my $checkPeptide = $origPeptide;
		$checkPeptide =~ s/[$modSymbols]//g;
		foreach my $site (keys %protSiteScores) {
			my ($pos) = $site =~ /[A-Z](\d+)/;
			if (substr($checkPeptide, $pos - 1, 2) eq "SP") {
				push (@SPsites, $site);
			}
		}
		if (@SPsites) {
			if (scalar(@SPsites) > $nSites) {
				## Check if any "SP" sites are overlapped with the original modification sites
				my $nOverlap = 0;
				foreach my $origSite (@origSites) {
					foreach my $SPsite (@SPsites) {
						if ($origSite eq $SPsite) {
							$nOverlap++;
						}
					}
				}
				if ($nOverlap > $nSites) {
					## If there are more SP sites (overlapped with the original modification sites) than $nSites,
					## we need to randomly select $nSites SP sites
					## Example,
					## Original peptide = SKLS#PS#PSLR
					## Candidate sites = S1, S4 and S6
					## PepSiteScores = (65.26, 65.26 and 63.80) for (S1, S4 and S6)
					## ProtSiteScores = (81.39, 80.40, 80.40) for (S1, S4 and S6)
					## - S1 is already selected based on protein-site scores
					## - One site need to be determined among S4 and S6
					## - Both of them are "SP" sites and overlapped with the original modification sites
					## - One of them needs to be "somehow" selected
					## - Then, go to the next outfile or peptide
					@SPsites = shuffle(@SPsites);
					for (my $i = 0; $i < scalar(@SPsites); $i++) {
						if ($i < $nSites) {
							$$outfileHash{'selectedSites'}{$SPsites[$i]} = $$outfileHash{'sites'}{$SPsites[$i]};
						} else {
							$$outfileHash{'alternativeSites'}{$SPsites[$i]} = $$outfileHash{'sites'}{$SPsites[$i]};
						}
					}
					next;
				} else {
					foreach my $origSite (@origSites) {
						for (my $i = 0; $i < scalar(@SPsites); $i++) {
							## Accept the "SP" site(s) matched to the original modification site(s)
							if ($origSite eq $SPsites[$i]) {
								$nSelected++;
								$$outfileHash{'selectedSites'}{$origSite} = $$outfileHash{'sites'}{$origSite};
								delete $protSiteScores{$origSite};
								splice (@SPsites, $i, 1);
							}
						}
					}
				}
				if ($nSelected == $nSites) {
					## Assign alternative sites
					## Then, go to the next outfile or peptide
					foreach my $site (keys %protSiteScores) {
						$$outfileHash{'alternativeSites'}{$site} = $$outfileHash{'sites'}{$site};
					}
					next;
				} else {
					## Among sites in the keys of %protSiteScores, ($nSites - $nSelected) sites need to be selected
					@SPsites = shuffle(@SPsites);
					for (my $i = 0; $i < scalar(@SPsites); $i++) {
						if ($i < $nSites - $nSelected) {
							$$outfileHash{'selectedSites'}{$SPsites[$i]} = $$outfileHash{'sites'}{$SPsites[$i]};
						} else {
							$$outfileHash{'alternativeSites'}{$SPsites[$i]} = $$outfileHash{'sites'}{$SPsites[$i]};
						}
					}
					$nRandom++;
					next;
				}
			} else {
				## Accept all "SP" sites
				foreach my $site (@SPsites) {
					$nSelected++;
					$$outfileHash{'selectedSites'}{$site} = $$outfileHash{'sites'}{$site};
					delete $protSiteScores{$site};
				}
				## Update the number of sites (to be determined)
				$nSites = $nSites - $nSelected;
				# Go to the next outfile if all sites have been determined
				## Otherwise, go to next step; look for and choose "S" sites
				next if ($nSites == 0);
			}
		}
		
		## Rule 2. Choose "S" site(s)
		$nSelected = 0;
		my @Ssites;
		foreach my $site (keys %protSiteScores) {
			if ($site =~ /^S/) {
				push (@Ssites, $site);
			}
		}
		if (@Ssites) {
			if (scalar(@Ssites) > $nSites) {
				## Check if any "SP" sites are overlapped with the original modification sites
				my $nOverlap = 0;
				foreach my $origSite (@origSites) {
					foreach my $Ssite (@Ssites) {
						if ($origSite eq $Ssite) {
							$nOverlap++;
						}
					}
				}
				if ($nOverlap > $nSites) {
					## If there are more S sites (overlapped with the original modification sites) than $nSites,
					## we need to randomly select $nSites SP sites
					## Example,
					## Original peptide = SRS#TS#EPEEAELSLSLAR
					## Candidate sites = S1, S3, T4 and S5
					## PepSiteScores = (49.91, 49.91, 49.90, 49.90) for (S1, S3, T4 and S5)
					## ProtSiteScores = (92.03, 49.99, 49.99, 49.99) for (S1, S3, T4 and S5)
					## - S1 is already selected based on protein-site scores
					## - One site need to be determined among S3 and S5
					## - Both of them are "S" sites and overlapped with the original modification sites
					## - One of them needs to be "somehow" selected
					## - Then, go to the next outfile or peptide
					@Ssites = shuffle(@Ssites);
					for (my $i = 0; $i < scalar(@Ssites); $i++) {
						if ($i < $nSites) {
							$$outfileHash{'selectedSites'}{$Ssites[$i]} = $$outfileHash{'sites'}{$Ssites[$i]};
						} else {
							$$outfileHash{'alternativeSites'}{$Ssites[$i]} = $$outfileHash{'sites'}{$Ssites[$i]};
						}
					}
					next;
				} else {
					foreach my $origSite (@origSites) {
						for (my $i = 0; $i < scalar(@Ssites); $i++) {
							## Accept the "S" site(s) matched to the original modification site(s)
							if ($origSite eq $Ssites[$i]) {
								$nSelected++;
								$$outfileHash{'selectedSites'}{$origSite} = $$outfileHash{'sites'}{$origSite};
								delete $protSiteScores{$origSite};
								splice (@Ssites, $i, 1);
							}
						}
					}
				}
				if ($nSelected == $nSites) {
					## Assign alternative sites
					## Then, go to the next outfile or peptide
					foreach my $site (keys %protSiteScores) {
						$$outfileHash{'alternativeSites'}{$site} = $$outfileHash{'sites'}{$site};
					}
					next;
				} else {
					## Among sites in the keys of %protSiteScores, ($nSites - $nSelected) sites need to be selected
					@Ssites = shuffle(@Ssites);
					for (my $i = 0; $i < scalar(@Ssites); $i++) {
						if ($i < $nSites - $nSelected) {
							$$outfileHash{'selectedSites'}{$Ssites[$i]} = $$outfileHash{'sites'}{$Ssites[$i]};
						} else {
							$$outfileHash{'alternativeSites'}{$Ssites[$i]} = $$outfileHash{'sites'}{$Ssites[$i]};
						}
					}
					$nRandom++;
					next;
				}
			} else {
				## Accept all "S" sites
				foreach my $site (@Ssites) {
					$nSelected++;
					$$outfileHash{'selectedSites'}{$site} = $$outfileHash{'sites'}{$site};
					delete $protSiteScores{$site};
				}
				## Update the number of sites (to be determined)
				$nSites = $nSites - $nSelected;
				# Go to the next outfile if all sites have been determined
				## Otherwise, go to next step; look for and choose "T" sites
				next if ($nSites == 0);
			}
		}
		
		## Rule 3. Choose "T" site(s)
		$nSelected = 0;
		my @Tsites;
		foreach my $site (keys %protSiteScores) {
			if ($site =~ /^T/) {
				push (@Tsites, $site);
			}
		}
		if (@Tsites) {
			if (scalar(@Tsites) > $nSites) {
				## Check if any "T" sites are overlapped with the original modification sites
				my $nOverlap = 0;
				foreach my $origSite (@origSites) {
					foreach my $Tsite (@Tsites) {
						if ($origSite eq $Tsite) {
							$nOverlap++;
						}	
					}
				}
				if ($nOverlap > $nSites) {
					@Tsites = shuffle(@Tsites);
					for (my $i = 0; $i < scalar(@Tsites); $i++) {
						if ($i < $nSites) {
							$$outfileHash{'selectedSites'}{$Tsites[$i]} = $$outfileHash{'sites'}{$Tsites[$i]};
						} else {
							$$outfileHash{'alternativeSites'}{$Tsites[$i]} = $$outfileHash{'sites'}{$Tsites[$i]};
						}
					}
					next;
				} else {
					foreach my $origSite (@origSites) {
						for (my $i = 0; $i < scalar(@Tsites); $i++) {
							## Accept the "T" site(s) matched to the original modification site(s)
							if ($origSite eq $Tsites[$i]) {
								$nSelected++;
								$$outfileHash{'selectedSites'}{$origSite} = $$outfileHash{'sites'}{$origSite};
								delete $protSiteScores{$origSite};
								splice (@Tsites, $i, 1);
							}
						}
					}
				}
				if ($nSelected == $nSites) {
					## Assign alternative sites
					## Then, go to the next outfile or peptide
					foreach my $site (keys %protSiteScores) {
						$$outfileHash{'alternativeSites'}{$site} = $$outfileHash{'sites'}{$site};
					}
					next;
				} else {
					## Among sites in the keys of %protSiteScores, ($nSites - $nSelected) sites need to be selected
					@Tsites = shuffle(@Tsites);
					for (my $i = 0; $i < scalar(@Tsites); $i++) {
						if ($i < $nSites - $nSelected) {
							$$outfileHash{'selectedSites'}{$Tsites[$i]} = $$outfileHash{'sites'}{$Tsites[$i]};
						} else {
							$$outfileHash{'alternativeSites'}{$Tsites[$i]} = $$outfileHash{'sites'}{$Tsites[$i]};
						}
					}
					$nRandom++;
					next;
				}
			} else {
				## Accept all "T" sites
				foreach my $site (@Tsites) {
					$nSelected++;
					$$outfileHash{'selectedSites'}{$site} = $$outfileHash{'sites'}{$site};
					delete $protSiteScores{$site};
				}
				## Update the number of sites (to be determined)
				$nSites = $nSites - $nSelected;
				# Go to the next outfile if all sites have been determined
				## Otherwise, go to next step; look for and choose "T" sites
				next if ($nSites == 0);
			}
		}
		
		## Rule 4. Choose "Y" site(s)
		$nSelected = 0;
		my @Ysites;
		foreach my $site (keys %protSiteScores) {
			if ($site =~ /^Y/) {
				push (@Ysites, $site);
			}
		}
		if (@Ysites) {
			if (scalar(@Ysites) > $nSites) {
				## Check if any "Y" sites are overlapped with the original modification sites
				my $nOverlap = 0;
				foreach my $origSite (@origSites) {
					foreach my $Ysite (@Ysites) {
						if ($origSite eq $Ysite) {
							$nOverlap++;
						}
					}
				}
				if ($nOverlap > $nSites) {
					@Ysites = shuffle(@Ysites);
					for (my $i = 0; $i < scalar(@Ysites); $i++) {
						if ($i < $nSites) {
							$$outfileHash{'selectedSites'}{$Ysites[$i]} = $$outfileHash{'sites'}{$Ysites[$i]};
						} else {
							$$outfileHash{'alternativeSites'}{$Ysites[$i]} = $$outfileHash{'sites'}{$Ysites[$i]};
						}
					}
					next;
				} else {
					foreach my $origSite (@origSites) {
						for (my $i = 0; $i < scalar(@Ysites); $i++) {
							## Accept the "Y" site(s) matched to the original modification site(s)
							if ($origSite eq $Ysites[$i]) {
								$nSelected++;
								$$outfileHash{'selectedSites'}{$origSite} = $$outfileHash{'sites'}{$origSite};
								delete $protSiteScores{$origSite};
								splice (@Ysites, $i, 1);
							}
						}
					}
				}
				if ($nSelected == $nSites) {
					## Assign alternative sites
					## Then, go to the next outfile or peptide
					foreach my $site (keys %protSiteScores) {
						$$outfileHash{'alternativeSites'}{$site} = $$outfileHash{'sites'}{$site};
					}
					next;
				} else {
					## Among sites in the keys of %protSiteScores, ($nSites - $nSelected) sites need to be selected
					@Ysites = shuffle(@Ysites);
					for (my $i = 0; $i < scalar(@Ysites); $i++) {
						if ($i < $nSites - $nSelected) {
							$$outfileHash{'selectedSites'}{$Ysites[$i]} = $$outfileHash{'sites'}{$Ysites[$i]};
						} else {
							$$outfileHash{'alternativeSites'}{$Ysites[$i]} = $$outfileHash{'sites'}{$Ysites[$i]};
						}
					}
					$nRandom++;
					next;
				}
			} else {
				## Accept all "Y" sites
				foreach my $site (@Ysites) {
					$nSelected++;
					$$outfileHash{'selectedSites'}{$site} = $$outfileHash{'sites'}{$site};
					delete $protSiteScores{$site};
				}
				## Update the number of sites (to be determined)
				$nSites = $nSites - $nSelected;
				if ($nSites > 0) {
					print "  Some sites still need to be determined in $peptide of $outfile\n";
					print LOGFILE "  Some sites still need to be determined in $peptide of $outfile\n";
					exit;
				}
			}
		}
	}
}
print "\n";
print LOGFILE "  Updated the modification-site scores in $nPeptides peptides (out of $nTotPeptides peptides)\n";

##########################################################
## Generate an output file (default = ID.lscore)	##
##########################################################
my $outputFile = $dir."/$Outfile";
open (OUT, ">", $outputFile);
open (IDmod, "<", $inputFile) or die "Cannot open $inputFile\n";
$header = <IDmod>;
print OUT $databaseHeader,"\n";
print OUT $header;

my $nOutputLines = 0;
while (<IDmod>) {
	$nOutputLines++;
	print "\r  Writing the result of $nOutputLines PSMs";
	chomp;
	my $line = $_;
	my @elems = split(/;/, $line);
	my $peptide = $elems[0];
	my $outfile = $elems[2];
	my ($AAstartPos) = $elems[14] =~ /AA(\d+)to/;
	
	my @origSiteInfoArray;
	for (my $i = 16; $i < scalar(@elems); $i = $i + 2) {
		my $origSiteInfo = $elems[$i].":".$elems[$i + 1];
		push (@origSiteInfoArray, $origSiteInfo);
	}
	my @siteInfoArray;
	my @sitePosArray;
	foreach my $site (keys %{$pepScores{$peptide}{$outfile}{'selectedSites'}}) {
		my ($siteAA, $sitePos) = $site =~ /([A-Z])(\d+)/; 
		my $siteScore = $pepScores{$peptide}{$outfile}{'selectedSites'}{$site};
		my $siteInfo = $siteAA.($sitePos + $AAstartPos - 1).":".$siteScore;
		push (@sitePosArray, $sitePos);
		push (@siteInfoArray, $siteInfo);
	}
	my $modPeptide = $peptide;
	$modPeptide =~ s/[$modSymbols]//g;
	my $nMethionines = 0;
	while ($modPeptide =~ /\@/g) {
		## n : the number of modification sites in the original peptide (=$nSites)
		$nMethionines++;
		my $pos = pos ($modPeptide) - $nMethionines - 2;
		push (@sitePosArray, $pos);
	}
	$modPeptide =~ s/\@//g;
	@sitePosArray = sort {$a <=> $b} @sitePosArray;
	my $nSites = 0;	
	for (my $i = 0; $i < scalar(@sitePosArray); $i++) {
		my $AA = substr($modPeptide, $sitePosArray[$i] - 1 + $nSites + 2, 1);
		if ($AA eq "S") {
			substr($modPeptide, $sitePosArray[$i] + $nSites + 2, 0, "#");
		} elsif ($AA eq "T") {
			substr($modPeptide, $sitePosArray[$i] + $nSites + 2, 0, "%");
		} elsif ($AA eq "Y") {
			substr($modPeptide, $sitePosArray[$i] + $nSites + 2, 0, "*");
		} elsif ($AA eq "M") {
			substr($modPeptide, $sitePosArray[$i] + $nSites + 2, 0, "@");
		}
		$nSites++;
	}
	my $outputLine = join(";", @elems[0..15]).";".join(",", @origSiteInfoArray).";".$modPeptide.";".join(",", @siteInfoArray);
	print OUT "$outputLine\n";
}
close (IDmod);
close (OUT);
print LOGFILE "  Wrote the result of $nOutputLines PSMs";