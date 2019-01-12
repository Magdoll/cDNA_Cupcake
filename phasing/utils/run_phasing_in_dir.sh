#!/bin/bash
DIR=$1
PLOIDY=2
MIN_BQ=13

echo "I'm going to $DIR"
STRAND=`grep "ref_strand=" $DIR/config|awk -F'=' '{print $2}'`
echo "strand is $STRAND"
PBID=`grep ">" $DIR/fake.fasta|awk -F"_" '{print $2}'`
PBID="fake_$PBID"

cat <<EOM >$DIR/run.sh

# 3. create mpileup

minimap2 -ax splice fake.fasta ccs.fastq > ccs.sam
#samtools view -b ../../all.fake.shortread.sorted.bam "$PBID" > mapped.shortread.bam

# 3. create mpileup
samtools view -bS ccs.sam > ccs.bam
samtools sort ccs.bam > ccs.sorted.bam
samtools mpileup --min-BQ $MIN_BQ -f fake.fasta -s ccs.sorted.bam > ccs.mpileup

#samtools merge out.bam mapped.shortread.bam ccs.sorted.bam
#samtools mpileup --min-BQ $MIN_BQ -f fake.fasta -s out.bam > ccs.mpileup

# 4. run phasing
run_phaser.py ccs.fastq ccs.sam ccs.mpileup fake.read_stat.txt fake.mapping.txt --strand $STRAND -o phased.nopartial -n $PLOIDY
run_phaser.py ccs.fastq ccs.sam ccs.mpileup fake.read_stat.txt fake.mapping.txt --partial_ok --strand $STRAND -o phased.partial -n $PLOIDY

# 5. create genome CCS BAM file for reference
#fq2fa.py ccs.fastq  # this can be removed in later 
#~/bin/gmap -D ~/share/gmap_db_new/ -d maize4 -f samse -n 0 -t 4 ccs.fasta > ccs.maize4.sam
#samtools view -bS ccs.maize4.sam > ccs.maize4.bam
#samtools sort ccs.maize4.bam > ccs.maize4.sorted.bam
#samtools index ccs.maize4.sorted.bam
EOM
