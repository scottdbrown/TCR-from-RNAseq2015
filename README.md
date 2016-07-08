Extracting TCR CDR3 sequences from RNAseq data.
===============================================

February 5, 2015
Scott Brown

Citation: Brown et al. Genome Medicine 2015, 7:125  doi:10.1186/s13073-015-0248-x
http://www.genomemedicine.com/content/7/1/125


--------------------------------------------------


The MiTCR tool does not use read pairing information, so all sequence data is first
pooled into a single fastq file. During creation of this pooled read file, reads
are checked that they only contain ACTGN bases, and are at least 40 bases in
length.

Space requirements:
This pipeline uses 10GB RAM.
This pipeline will create a copy of the sequence data, so will double the space requirement.

Pipeline:

## All fastq files for a single sample, unzipped and uncompressed, must reside in
## a single directory.

~$ python2 combineAndCleanFastq.py /path/to/fastq/files/ /working/directory/

## This will output /working/directory/reads.fq and /working/directory/numReads.txt
## reads.fq: Cleaned, combined RNA-seq file for input into MiTCR.
## numReads.txt: The number of RNA-seq reads to be used for MiTCR.

## To run MiTCR on RNA-seq data, the correct parameter set must be selected. This
## is based on TCR chain and read length. Parameter files are in the format:
## rnaseq[READ_LEN]bp[TCR_CHAIN].pset
## Choose READ_LEN (50, 76, 101) that is closest to actual read length
## Choose TCR_CHAIN (TRA, TRB) that matches chain you are interested in.

~$ java -Xmx10g -jar mitcr.jar -pset rnaseq50bpTRA.pset /working/directory/reads.fq /working/directory/outputTRA50bp.txt

## And the other chain...

~$ java -Xmx10g -jar mitcr.jar -pset rnaseq50bpTRB.pset /working/directory/reads.fq /working/directory/outputTRB50bp.txt

## Output files contain all CDR3 sequences that were extracted.

