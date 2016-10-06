Extracting TCR CDR3 sequences from RNAseq data.
===============================================

February 5, 2015
Scott Brown

Citation: Brown et al. Genome Medicine 2015, 7:125  doi:10.1186/s13073-015-0248-x
http://www.genomemedicine.com/content/7/1/125


--------------------------------------------------

*NOTE:* See [below](#simulateTCR) for information on simulating TCR transcripts.

The MiTCR tool does not use read pairing information, so all sequence data is first
pooled into a single fastq file. During creation of this pooled read file, reads
are checked that they only contain ACTGN bases, and are at least 40 bases in
length.

Space requirements:
This pipeline uses 10GB RAM.
This pipeline will create a copy of the sequence data, so will double the space requirement.

Pipeline:

All fastq files for a single sample, unzipped and uncompressed, must reside in a single directory.

`~$ python2 combineAndCleanFastq.py /path/to/fastq/files/ /working/directory/`

This will output /working/directory/reads.fq and /working/directory/numReads.txt
reads.fq: Cleaned, combined RNA-seq file for input into MiTCR.
numReads.txt: The number of RNA-seq reads to be used for MiTCR.

To run MiTCR on RNA-seq data, the correct parameter set must be selected. This
is based on TCR chain and read length. Parameter files are in the format:
rnaseq[READ_LEN]bp[TCR_CHAIN].pset

Choose READ_LEN (50, 76, 101) that is closest to actual read length
Choose TCR_CHAIN (TRA, TRB) that matches chain you are interested in.

*NOTE*: mitcr.jar obtained from [github.com/milaboratory/mitcr/releases](https://github.com/milaboratory/mitcr/releases)

`~$ java -Xmx10g -jar mitcr.jar -pset rnaseq50bpTRA.pset /working/directory/reads.fq /working/directory/outputTRA50bp.txt`

And the other chain...

`~$ java -Xmx10g -jar mitcr.jar -pset rnaseq50bpTRB.pset /working/directory/reads.fq /working/directory/outputTRB50bp.txt`

Output files contain all CDR3 sequences that were extracted.


<a name="simulateTCR"></a>
###Simulating TCR transcripts.

File: [simulate_TCR_transcripts.py](simulate_TCR_transcripts.py)

Requires Python 2.

Usage:
`python2 simulate_TCR_transcripts.py /path/to/references/ outputFile.fq numberOfAlphaBetaPairs`

Directory containing references expects a separate file for each gene, ex TRAV.fa, TRAJ.fa, TRAC.fa, etc. Different alleles for a gene are contained within the file, and have the entire sequence on a single line.

Quality scores for fastq are set uniformly high ("J").

TCR genes are selected randomly from a uniform distribution. Junctional nucleotide additions and subtractions are selected based on observed frequencies.
