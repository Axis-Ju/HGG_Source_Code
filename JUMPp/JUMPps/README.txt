----------------------
Contents of this file
----------------------

 * Introduction
 * Software Requirements
 * Hardware Requirements
 * Installation
 * Command Line Arguments
 * Maintainers

----------------------
Introduction
----------------------

JUMP is a new hybrid database search algorithm that generates amino acid tags and ranks peptide spectrum matches by the tags and pattern matching. JUMP can use all potential sequence tags, as short as only one amino acid, in database search. Thus JUMP has high sensitivity and specificity to differentiate true matches from false matches. The program also provides additional features, such as (i) identification of co-eluted peptides from mixture MS/MS spectra, and (ii) assignment of modification sites by tags. The current version is optimized for the analysis of high resolution MS/MS spectra. The program is written in perl and is designed for high performance parallel computing systems.
 
----------------------
Software Requirements
---------------------- 

The program is written in Perl. It should run on Linux system with a Perl5. The minimum required Perl version should be Perl 5.8 or better.

Essentially, It can be executed without installing any Perl dependencies.
 
----------------------
Hardware Requirements
---------------------- 

The program can be run on either high performance computing systme or a single server. 
 
To run on a cluster:
 Batch-queuing system: SGE, version 6.1u5, qsub and qstat
 32 GB memory on each node

To run on a single server
  32 GB memory
  2 GHz CPU processors with a minimum of 4 cores
  
----------------------
Installation
---------------------- 

After downloading the source code, you can put it in any working directory (e.g. /home/xxxx/JUMP). 
**If you install the program on the cluster, the folder containing all source code is required to be accessible to each node. For example, you can put it on your local home directory. 


----------------------
Command Line Arguments
----------------------

Example:  perl /home/xxxx/JUMP/jump.pl -p <JUMP parameter file> <MS/MS data file(s) with .mzXML>
-p <file> specifies a JUMP parameter file
<MS/MS data file(s)> specifies one or multiple MS/MS data file(s)

The mzXML file is converted from .RAW file by various programs, such as ReAdW or msconvert. When using msconvert, please select 32-bit for the option of "Binary encoding precision" and unselect the option of "use zlib compression".

An example parameter file is provided in the package of source code. 

----------------------
Maintainers
----------------------

JUMPps was published: https://www.stjuderesearch.org/site/lab/peng

Please contact Xusheng Wang (xusheng.wang@stjude.org) and Junmin Peng (junmin.peng@stjude.org) for software support.