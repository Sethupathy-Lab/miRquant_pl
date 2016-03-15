# set path
#
export PATH=/proj/.test/roach/miRNA/bin:./:$PATH
export PYTHONPATH=/proj/.test/roach/miRNA/lib/python/
export PATH=$PATH:/proj/.test/roach/miRNA/SHRiMP_2_2_2/bin
export SHRIMP_FOLDER=/proj/.test/roach/miRNA/SHRiMP_2_2_2

export JB_GENOME_DIR=/proj/seth_lab/projects/genome/


module unload bedtools
module unload perl 
module unload git

module load bedtools/2.17.0
module load ghostscript/8.71
module load r/2.15.1
module load git/1.8.5.3
module load star/2.3.0
module load muscle/3.8.31
module load bowtie/1.1.0
module load bowtie2/2.2.1
module load fastqc/0.11.2
module load gcc/4.8.1
module load gnuplot/4.4.0
module load samtools/1.3
module load bedtools/2.17.0
module load mathematica/9.0
module load picard/1.88
module load bamtools/1.0.2
module load perl/5.12.0
module load tophat/2.0.11
module load boost/1.53.0
