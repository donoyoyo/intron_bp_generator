
Perl script to generate sequence information regarding splice sites and branch points.
It was designed for yeast, so uses a strict definition of branchpoint

searching in this order

1.) ACTAA

2.) [AG][TC]T[AG]A[TC][AG]
aka PYTPAYP

3.) [TC]T[AG]A[TC][ACTG]
aka YTPAY


USAGE:

./intron_branchpoint_gen.pl inputintron.bed annotation.bed genome.2bit outputfile.tsv conservation.bw

inputintron.bed contains bed format from data that has introns
annotation.bed contains the annotated genes in the genome
genome.2bit is the genome sequence in 2bit format
outputfile is a tab seperated value text file sorted by total bed score
conservationfile (bigwig format , optional)

You will need twoBitToFa in your path

(and hgBigWigSummary if you are using conservation)

useful files
ares_sc3_isoforms.bed
intron_branchpoint_gen.pl
sacCer3.2bit
sacCer3.phastCons7way.bw

using those files:
 
 ./intron_branchpoint_gen.pl introndata.bed ares_sc3_isoforms.bed sacCer3.2bit outputfile.tsv sacCer3.phastCons7way.bw

The phastCons big wig file is available several places including here:

http://hgdownload.cse.ucsc.edu/goldenPath/sacCer3/phastCons7way/sacCer3.phastCons7way.bw


 sacMik2
 sacBay2
 ares_sc3_isoforms.bed
 sc3intronsfeb15.bed


