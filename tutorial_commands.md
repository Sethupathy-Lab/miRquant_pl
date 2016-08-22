## miRquant tutorial commands

#### First things first
The other readme is more expansive with the decriptions of the steps and explanations of the commands.  
However, this will allow you to copy and paste the different steps and complete the tutorial, then try the tutorial using the main directions.

Remember:  
Tab completion is the most useful of all skills to learn.  It'll speed up your navigation of the command line and avoid spelling issues.  
The * is a wildcard, meaning * alone is all files, m* is all files starting with m, *txt are all files ending with txt, and so forth.

#### Set path to tutorial folder to a variable
This section will be the only section different for you for the rest of the tutorial.
```
tutPath='/proj/seth_lab/users/ONYEN/miRquant_tutorial'

where ONYEN is your ONYEN
```

The tutorial path is the path to the location of the miRquant tutorial folder, which contains miRquant and tut_samples.
Let's make sure we have that correct by typing:
```
echo $tutPath
```
Then change into the miRquant tutorial folder.
```
cd $tutPath/miRquant
```

#### Load environmental variables
From here down, you should be able to copy and paste commands.
```
source uncENV.sh

Output should look like this:

	  Genome indexes for use with BOWTIE are available
	    in /proj/seq/data .

	  bwa/bowtie/bowtie2 indexes are located in "Sequence" directory
	    under the top level genome directory with name in CAPS, e.g. MM9_UCSC
```

#### Generate the adapter files
Change into the scripts folder and run the generate_adapter_files.py
```
cd scripts

python generate_adapter_files.py $tutPath/tut_samples/*.fastq

Output should look like this:

Sample: SampleA_ATCACG_tut.fastq
Output saved as: /proj/seth_lab/users/Matt/miRquant_tutorial/tut_samples/SampleA_ATCACG_tut.adaptor
Adaptor sequence: TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG
```

#### Run chainSubmission.sh
Move up one directory to the location of the chainSubmission.sh file.
```
cd ..
```
Run chainSubmission.sh
```
bsub -oo chainSubmission.log bash chainSubmission.sh 10 1 hsa 33 NoGS $tutPath/tut_samples/SampleA_ATCACG_tut.fastq
```
The test sample is from human, and that is why is hsa given (which could be mmu for mouse or rno for rat).

Wait for all your jobs to finish running before continueing (you can check this by typing bjobs)

#### Run runC.sh

You might want to check out the chainSubmission.log see how the run went.  You can do this by typing:
```
more chainSubmission.log

space key will cycle you down
q will exit you out
```

Run runC.sh
```
bash runC.sh hsa $tutPath/tut_samples/*/IntermediateFiles/g1Results/CHR*.results
```

#### Run post_runC.sh
This script concatenates the results of runC.sh into single files.  To run it, type:
```
bash post_runC.sh $tutPath/tut_samples/*/IntermediateFiles/g1Results/
```

#### Run process_all_summary2tab.pl
Change to the scripts folder
```
cd scripts
```
and run process_all_summary2tab.pl by typing:
```
bsub -oo summary2tab.log perl process_all_summary2tab.pl $tutPath/miRquant hsa $tutPath/tut_samples/*/IntermediateFiles/g1Results/shift_summary.txt
```

## Final processing
#### Generate mapping table
```
python generate_mapping_info.py $tutPath/tut_samples/*/*.stats
```
The output file is called MappingInfoTable.csv and looks like this:
```
cat MappingInfoTable.csv

Sample_name,SampleA_ATCACG_tut.
File,/proj/seth_lab/users/Matt/miRquant_tutorial/tut_samples/SampleA_ATCACG_tut./SampleA_ATCACG_tut.stats
Total Reads,100000
Trimmed Reads,90252
Percent Trimmed Reads,90.25
Too Short Reads,7012
Percent Too Short,7.01
Exact Match Reads,58610
Percent Exact Matches,64.94
Mismatch Reads,31642
Percent Mismatched,35.06
Mapped Reads, 67169
Percent Mapped,74.42
miR Mapped Reads, 45760
Percent miR Mapped,68.13
tRNA Mapped Reads, 527
Percent tRNA Mapped,0.78
```

#### Generate length histogram
```
python lenDist.py $tutPath/tut_samples/*/IntermediateFiles/*_O10_E1.fq --image
```
The outputs are lenDist.tsv and lenDistHistogram.png (image)

lenDist.tsv should look like this:
```
cat lenDist.tsv

	SampleA_ATCACG_tut._O10_E1.fq
14	0.0292957496787
15	0.0356778797146
16	0.0226698577317
17	0.0270686522182
18	0.0372512520498
19	0.0352014359793
20	0.039500509684
21	0.123565128751
22	0.264758675708
23	0.135996986216
24	0.0442981872978
25	0.0204095200106
26	0.0180494615078
27	0.016609050215
28	0.0175508575987
29	0.0131853033728
30	0.0151021584009
31	0.0216504897398
32	0.0128418206799
33	0.014758675708
34	0.0118556929486
35	0.010902805478
36	0.00699153481363
37	0.00591676638745
38	0.00717989629039
39	0.0055178832602
40	0.00278110180384
41	0.00341266675531
```

#### Generate RPM files
```
python genNormalRPM.py hsa $tutPath/tut_samples/*/TAB_lenDist_summary.txt
```

The output of this will be RPM_all.csv, RPM_miRs_only.csv, RPM_miRs_over_100.csv.  
The top ten lines of these files look similar to this:
```
head RPM*

==> RPM_all.csv <==
,SampleA_ATCACG_tut.
CHR1:228750515:M ,467.476949646
CHR1:228777350:M ,200.626497282
CHR3:32027818:M ,187.586291896
CHR1:228754993:M ,299.244798977
CHR1:228757210:M ,456.92014534
CHR12:127650532:M ,643.472025093
CHR1:228754991:M ,583.601797011
CHR17:33478293:M ,163.765810386
CHR17:33478295:M ,506.185232101

==> RPM_miRs_only.csv <==
,SampleA_ATCACG_tut.
hsa-mir-423-3p,163.765810386
hsa-mir-181a-2-5p,357.307222659
hsa-mir-26a-1-5p,7974.89870575
hsa-mir-192-5p,35834.9368726
hsa-mir-203,976.026714838
hsa-mir-199b-3p,506.185232101
hsa-mir-425-5p,208.429213218
hsa-mir-25-3p,1161.24847364
hsa-mir-31-5p,1012.3704642

==> RPM_miRs_over_100.csv <==
,SampleA_ATCACG_tut.
hsa-mir-423-3p,163.765810386
hsa-mir-181a-2-5p,357.307222659
hsa-mir-26a-1-5p,7974.89870575
hsa-mir-192-5p,35834.9368726
hsa-mir-203,976.026714838
hsa-mir-199b-3p,506.185232101
hsa-mir-425-5p,208.429213218
hsa-mir-25-3p,1161.24847364
hsa-mir-31-5p,1012.3704642
```

#### Generate the RPMMM file
```
python genNormalRPMMM.py hsa $tutPath/tut_samples/*/TAB_3p_summary_miR.txt
```

The output from this should be RPMMM.csv.
