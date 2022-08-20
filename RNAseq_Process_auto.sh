#!/bin/bash

help(){
	cat << EOF
Description:
	This mini tool can be used to process clean RNA-seq data and produce bam, bigwig, bedgraph files.
	
	CAUTION: (1) THE INPUT DATA MUST BE CLEANED!!!
		 (2) Until now, this tool can only deal with strand specific RNA-seq data (ssRNA-seq).
Required tools:
	- STAR
	- bedtools
	- samtools
	- featureCounts
Usage:
	RNAseq_Process_auto.sh -1 <r1.fq> -2 <r2.fq> -r <ref> -g <gtf> -b <bed> -n <processers> -p <prefix> -o <output>
Options:
	-1	the path of read1 fastq, can't be gzip. [PATH]
	-2	the path of read2 fastq, can't be gzip. [PATH]
	-r	the ref genome for STAR to do reads alignments. [PATH]
	-g	the GTF annotation file for reads summary. [PATH]
	-b	the genome annotation bed file used for determine library type. [PATH]
	-n	the number of processers to perform this task, at least 2. [NUM]
	-p	the prefix of output file. [CHAR]
	-o	the output dir for this task, can be created automatically. [CHAR]
Owner:
	Kunming Shui, skm@smail.nju.edu.cn, Nanjing University.
Last Modified:
	2022-07-20
Modified Log:
	2022-07-18	RNAseq_Process_auto.sh first released.[v 0.1.0]
	2022-07-20	Modified a bug about the strand information determination.[v 0.1.1]
EOF
}

#version information
version="v 0.1.1"

if [ $# -eq 0 ] || [ $1 = "-h" ] || [ $1 = "--help" ]
then
	help
	exit 1
fi

if [ $1 = "-v" ] || [ $1 = "--version" ]
then
	echo $version
	exit 1
fi

echo "=================== Begin Analysis at $(date) ==================="

while getopts "1:2:r:g:b:n:p:o:" option
do
	case $option in
		1) read1=$OPTARG;;
		2) read2=$OPTARG;;
		r) ref=$OPTARG;;
		g) gtf=$OPTARG;;
		b) bed=$OPTARG;;
		n) num=$OPTARG;;
		p) prefix=$OPTARG;;
		o) out=$OPTARG;;
	esac
done

echo "=================== Reads Alignment with STAR ==================="

if [ -d $out ]
then
	echo "$out is exiting."
else
	mkdir -p $out
fi

#data check
echo "The read1 file is $read1"
echo "The read2 file is $read2"

STAR --runThreadN $num \
	--genomeDir $ref \
	--readFilesIn $read1 $read2 \
	--outSAMtype BAM SortedByCoordinate \
	--outBAMsortingThreadN $[ $num/2 ] \
	--outFileNamePrefix $out/$prefix

echo "=================== Filter out low quality alignments ==================="
samtools view -q 255 -b -@ $num -o $out/$prefix.highquality.bam $out/${prefix}Ali*.bam
echo "Finished..."

echo "=================== Determine the library type ==================="
class=`infer_experiment.py -i $out/$prefix.highquality.bam -r $bed`
num1=`echo "$class" | grep 'Fraction of reads explained by "1++,1--,2+-,2-+":' | cut -d " " -f7`
num2=`echo "$class" | grep 'Fraction of reads explained by "1+-,1-+,2++,2--":' | cut -d " " -f7`

if [ `echo "$num1 > $num2" | bc` -eq 1 ]
then
	library="S"
else
	library="R"
fi

if [ $library = "S" ]
then
	echo "Your RNA-seq library is stranded..."
else
	echo "Your RNA-seq library is reverse stranded..."
fi

cd $out

samtools index $prefix.highquality.bam
samtools view -f 99 -b -@ $num -o $prefix.r1.fwd.bam $prefix.highquality.bam
samtools view -f 163 -b -@ $num -o $prefix.r2.fwd.bam $prefix.highquality.bam
samtools view -f 147 -b -@ $num -o $prefix.r2.rev.bam $prefix.highquality.bam
samtools view -f 83 -b -@ $num -o $prefix.r1.rev.bam $prefix.highquality.bam

if [ $library = "S" ]
then
	samtools merge -@ $num tmp.forward.bam $prefix.r1.fwd.bam $prefix.r2.rev.bam
	samtools merge -@ $num tmp.reverse.bam $prefix.r1.rev.bam $prefix.r2.fwd.bam
else
	samtools merge -@ $num tmp.forward.bam $prefix.r1.rev.bam $prefix.r2.fwd.bam
	samtools merge -@ $num tmp.reverse.bam $prefix.r1.fwd.bam $prefix.r2.rev.bam
fi

echo "=================== Perform reads summary with featureCounts ==================="
if [ $library = "S" ]
then
	type=1
else
	type=2
fi

mkdir featureCounts
featureCounts -a $gtf \
	-o featureCounts/${prefix}.txt \
	-t exon \
	-g gene_id \
	-s $type \
	-p \
	--countReadPairs \
	-B \
	-C \
	-T $num  $prefix.highquality.bam

echo "=================== Produce bedgraph, bigwig files ==================="
samtools sort -@ $num -o $prefix.fwd.bam tmp.forward.bam
samtools sort -@ $num -o $prefix.rev.bam tmp.reverse.bam

rm tmp.forward.bam tmp.reverse.bam
samtools index $prefix.fwd.bam
samtools index $prefix.rev.bam

#bedgraph file
genomeCoverageBed -bga -ibam $prefix.fwd.bam -split > $prefix.fwd.bedgraph
genomeCoverageBed -bga -ibam $prefix.rev.bam -split > $prefix.rev.bedgraph

#bigwig file
bamCoverage -b $prefix.fwd.bam -o $prefix.fwd.bigwig -bs 10 --normalizeUsing CPM -p $num
bamCoverage -b $prefix.rev.bam -o $prefix.rev.bigwig -bs 10 --normalizeUsing CPM -p $num

echo "$class" >> Process.log
echo $library >> Process.log

cd -

echo "=================== Finished at $(date) ==================="
