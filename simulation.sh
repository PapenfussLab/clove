#!/bin/bash


java -cp bin VariationRandomization chr21_22.fa.fai 500000 human_test hom

java -cp lib/commons-lang3-3.1.jar:lib/picard-1.77.jar:lib/sam-1.77.jar:bin ReferenceSimulation chr21_22.fa human_test.fa

java  -jar -Xmx2048m ~/SimSeq/SimSeqNBProject/store/SimSeq.jar -1 100 -2 100 --error ~/SimSeq/examples/hiseq_mito_default_bwa_mapping_mq10_1.txt --error2 ~/SimSeq/examples/hiseq_mito_default_bwa_mapping_mq10_2.txt -l 300 -s 30 -n 10500000 -r human_test.fa_reference.fa -o _human_test_reads.sam -u 0.01

SamToFastq I=_human_test_reads.sam FASTQ=human_test_reads_R1.fastq F2=human_test_reads_R2.fastq VALIDATION_STRINGENCY=SILENT

bowtie2 --local -p 16 -x chr21_22_bt2 -1 human_test_reads_R1.fastq -2 human_test_reads_R2.fastq | samtools view -bS -o human_test_reads.bam -

samtools sort human_test_reads.bam human_test_reads_s

samtools index human_test_reads_s.bam


