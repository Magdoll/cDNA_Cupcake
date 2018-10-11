#!/bin/bash

#source /mnt/software/Modules/current/init/bash
#module load blast/2.6.0+


INPUT=$1
#DB=/pbi/dept/secondary/siv/gconcepcion/db/ncbi/refseq/complete
DB=/pbi/dept/secondary/siv/gconcepcion/db/ncbi/nt
OUTPUT=${INPUT}.blastn


blastn -task megablast -db ${DB} -query ${INPUT} -evalue 0.1 -num_alignments 1 -parse_deflines -outfmt '7 qaccver saccver stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore' -num_threads 30 >${OUTPUT}
