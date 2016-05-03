---
...

**Small RNA Pipeline Setup 05/03/16**

\

**Check out a copy of the smRNA pipeline code:**

If you don’t already have a directory:

\

\$ mkdir /proj/seth\_lab/users/ONYEN

\$ cd /proj/seth\_lab/users/ONYEN

\$ module load git

\$ git clone https://github.com/Sethupathy-Lab/smRNA\_pipeline.git

\

Due to size constraints, the resources folder couldn't be hosted on
github. These can be generated using the scripts and commands in the
resources folder, or pre-generated resource files can be obtained by
emailing [Matt.Kanke@gmail.com](mailto:Matt.Kanke@gmail.com)

\

Now you will have a directory:
/proj/seth\_lab/users/ONYEN/smRNA\_pipeline

Most code is run from this directory!

\

\

**Make a folder for your run in the smallRNA directory**
(/proj/seth\_lab/projects/smallRNA/)

\

\$ mkdir /proj/seth\_lab/projects/smallRNA/MY\_PROJECT\_NAME

\

**Copy files from where ever they are to your smallRNA directory:**

\$ cp /path/to/sequencing\_files/\*
/proj/seth\_lab/projects/smallRNA/MY\_PROJECT\_NAME

\

**Change to project directory and uncompress files:**

\

\$ cd /proj/seth\_lab/projects/smallRNA/MY\_PROJECT\_NAME

\

\$ bsub gunzip \*.gz

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

\

**Make the adaptor files for the samples:**

\$ cd /proj/seth\_lab/users/ONYEN/smRNA\_pipeline/scripts

\$ module load python/2.7.6

\$ bsub -oo adapter.log python generate\_adapter\_files.py

/proj/seth\_lab/projects/smallRNA/MY\_PROJECT\_NAME/\*.fastq

\

For file names such as File\_Cond\_AAGGCC\_etc.fastq: the following
script will pull the barcode out of the file name and create the adapter
files. It should automatically determine where the barcode is in the
file name, but if it isn't able to (which will be indicated by an error
in the log), you need to tell the script where the index is. It
separates the parts joined by the “\_” character. In the above example
the word 'File' is in position 0, the word 'Cond' is in position 1, and
the index AAGGCC is in position 2. So the command would be:

\

\$ bsub -oo adapter.log python generate\_adapter\_files.py 2
/proj/seth\_lab/projects/smallRNA/MY\_PROJECT\_NAME/\*.fastq

\

\

**Double check your adapter files have generated correctly:**

\

\$ cat adapter.log

\

For each sample, you should see something like this:

\

Sample: SAMPLE.fastq

Output saved as:
/proj/seth\_lab/projects/smallRNA/MY\_PROJECT\_NAME/SAMPLE.adaptor

Adaptor sequence:
TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACGATCAGATCTCGTATGCCGTCTTCTGCTTG

\

Your adapter sequence should only include nucleotides, if it doesn't,
try re-running the script and defining where the barcode is in the
sample name. The adapter files are saved in the same location as the
sample fastq.

\

ALTERNATLY: If not a TrueSeq library prep you need a file with the same
name as the .fastq file with the extension .adaptor (note the spelling!
There was originally two separate files .adapter and .adaptor, and I
have never bothered to change the code.): The file must contain one line
with the 3’ adapter for the prep used.

\

\

**Load proper modules and environmental variables:**

\$ cd /proj/seth\_lab/users/ONYEN/smRNA\_pipeline

\$ module clear

Answer yes to the prompt on whether to clear modules, then:

\$ source uncENV.sh

\

\

\

\

\

\

**Run the chain submission:**

From the smallRNApipeline directory:

\$ cd /proj/seth\_lab/users/ONYEN/smRNA\_pipeline

\

Then:

\$ bsub -o chainSub.log bash chainSubmission.sh 10 1 hsa 33 NoGS
/proj/seth\_lab/projects/smallRNA/MY\_PROJECT\_NAME/\*.fastq

\

chainSub.log – log file can be named anything, but chainSub.log makes
sense

hsa – hsa for human; mmu for mouse, rno for rat

\

Review /proj/seth\_lab/users/ONYEN/smRNA\_pipeline/example.README for an
explanation of the inputs.

\

Check that your jobs are running

\$ bjobs

\

Once all jobs have finished ( “No unfinished jobs” message):

This might take awhile depending on the number of samples run.

Once the chain submission has finished, you can see if there were any
errors in your logfile (chainSub.log in the command above).

\

In your MY\_PROJECT\_NAME directory, there will be a directory for each
FILENAME.fastq (called FILENAME.)

\

In that directory, there will be an IntermediateFiles subdirectory and a
FILENAME.stats file.

\

To get an idea of how the trimming and aligning is looking, type:

\

\$ cat /proj/seth\_lab/projects/smallRNA/MY\_PROJECT\_NAME/\*/\*.stats

\

file:/proj/seth\_lab/projects/smallRNA/MY\_PROJECT\_NAME/FILENAME.fastq

TotReads:21189608.00000000000000000000

TrimmReads:17621310.00000000000000000000

ShortReads:2500007.00000000000000000000

EMhits:10710591

EMmiss:6910719

\

TotReads = Total number of reads for this file

TrimmReads = \# of reads successfully trimmed of 3’ adapter

ShortReads = \# of reads too short after trimming (\< 13 )

EMhits = \# of reads with an exact alignment to the genome

EMmiss = \# of reads that fail to exactly align to genome

\

\

\

\

**Run the next stage to collect results:**

From your pipeline directory
(/proj/seth\_lab/users/ONYEN/smRNA\_pipeline):

\$ bash runC.sh hsa
/proj/seth\_lab/projects/smallRNA/MY\_PROJECT\_NAME/\*/IntermediateFiles/g1Results/CHR\*.results

\

Once all of those jobs have finished running, run:

\$ bash post\_runC.sh
/proj/seth\_lab/projects/smallRNA/MY\_PROJECT\_NAME/\*/IntermediateResults/g1Results

\

**Run the next stage to generate TAB separated files:**

\$ cd /proj/seth\_lab/users/ONYEN/smRNA\_pipeline/scripts

\$ bsub –o logFileName.log perl process\_all\_summary2tab.pl
/proj/seth\_lab/users/ONYEN/smRNA\_pipeline hsa
/proj/seth\_lab/projects/smallRNA/MY\_PROJECT\_NAME/\*/IntermediateFiles/g1Results/shift\_summary.txt

\

After run finishes, you should see:

\$ cd /proj/seth\_lab/projects/smallRNA/MY\_PROJECT\_NAME/

\$ cat \*/\*.stats

\

file:/proj/seth\_lab/projects/smallRNA/MY\_PROJECT\_NAME/FILE.fastq

TotReads:6149484.00000000000000000000

TrimmReads:3730081.00000000000000000000

ShortReads:1938220.00000000000000000000

EMhits:2509867

EMmiss:1220214

Mapped: 2881057.09082651

miRMapped: 1104952.30513136

\

Mapped and miRMapped indicate the number of reads mapped to the genome
and to miRNAs respectively.

\

The log file above will also contain a table that you can use to put
together the mapping stats for your project.

\

**TCGA files have adapters already trimmed.**

There is no need to generate the adapter files in the previous section.

All following steps remain the same, but use the chainSubmissionTCGA
script.

\

**Outputs:**

For each Sample:

TAB\_3p\_summary.txt - 3'-end differences

TAB\_3p\_summary\_miR.txt - 3'-end differences (miRNA loci only)

TAB\_ed\_summary.txt - central differences

TAB\_lenDist\_summary.txt - length differences

Shrimp\_results.bed - bed file containing all results

\

\

\

\

\

**To generate****mapping summary file****:**

\$ cd /proj/seth\_lab/users/ONYEN/smRNA\_pipeline/scripts

\$ bsub –o map\_summ.log perl generate\_mapping\_info\_table.pl
relative\_path\_to/MY\_PROJECT\_NAME/\*/\*.stats

\

Output:

MappingInfoTable.tsv - Mapping summary for all samples

\

We shoot for short reads around 10%, Trimmed reads around 85%, mapped at
70-80%, miR mapped around 70%.

\

**To generate length distribution file:**

\$ cd /proj/seth\_lab/users/ONYEN/smRNA\_pipeline/scripts

\$ module load python/2.7.6

\$ module load r

\$ bsub –o lenDist.log python lenDist.py
/proj/seth\_lab/projects/smallRNA/MY\_PROJECT\_NAME/\*/IntermediateFiles/\*O10\_E1.fq
--image

\

Output:

lenDist.tsv - length distribution for the fastq, should see a bump at
\~22 for

miRNAs, and maybe a bump around 30 for tRNAs

lenDistHistogram.png - *optional*, length distribution image, output if
--image flag included

\

**To generate normalized expression (RPMM) across multiple samples:**

\$ cd /proj/seth\_lab/users/ONYEN/smRNA\_pipeline/scripts

\$ module load python/2.7.6

\$ bsub –o RPMM.log python genNormalRPM.py hsa
/proj/seth\_lab/projects/smallRNA/MY\_PROJECT\_NAME/\*/TAB\_lenDist\_summary.txt

\

Outputs:

RPM\_all.tsv - RPM for everything

RPM\_miRs\_only.tsv - RPM for miRs

RPM\_mirRs\_over\_100.tsv - RPM for miRs if the RPM for the miR was over
100 for at

least sample

\

**To****generate sample correlations****:**

Correlation between samples in the miRs over 100 output. This should be
possible to do in excel. I have included the R script I use, which
requires variables to be switched within and is run on my local machine.
The R script is smRNAseq\_correlation.R

\

Output from R script:

sample\_correlation\_values.tsv - Table of correlations between samples

Sample\_correlation\_heatmap.png - Heatmap image visualizing data in
correlation table

\

\

The above steps create single sheets to be assembled for the final excel
report.

\

The final tab is done manually in excel to get the fold-change and
p-value between sample groups.
