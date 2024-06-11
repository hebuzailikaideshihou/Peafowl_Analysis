#!/bin/bash

BuildDatabase -name Pavonew  -engine ncbi White_Pavo.fa

RepeatModeler -engine ncbi -database Pavonew  -threads 56 -LTRStruct

RepeatMasker -pa 32 -poly -html -gff -a -lib all_White_Pavo_final.fasta -dir ./ White_Pavo.fa


braker.pl --genome=White_Pavo.fa.masked --prot_seq=all.fa --threads=8 --gff3 --min_contig=4000

gth -genomic genome.fa -protein all.fasta -gff3out -intermediate -o genomethreader.out
perl ./releated_scripts/Convert_gth.pl genomethreader.out > genomethreader.gff


hisat2-build -p 28 genome.fa genome &>bowtie2-build_run.log
cut -d \",\" -f 1 $RNA_Seq |sed 's/$/.gtf/' >mergelist.txt #prepare a list for stringtie-merge


hisat2 -p 28 -x genome --new-summary -S White_PavoCA-Y1_RRAS22634-V.sam -1 /gpfs/home/yexinhai/project/16_White_Pavo_genome/01_data/02_rna_seq/White_PavoCA-Y1_RRAS22634-V_1.clean.fq -2 /gpfs/home/yexinhai/project/16_White_Pavo_genome/01_data/02_rna_seq/White_PavoCA-Y1_RRAS22634-V_2.clean.fq
samtools sort -@ 8 -o White_PavoCA-Y1_RRAS22634-V.bam White_PavoCA-Y1_RRAS22634-V.sam
stringtie -p 28 -o White_PavoCA-Y1_RRAS22634-V.gtf White_PavoCA-Y1_RRAS22634-V.bam

stringtie --merge -p 28 -o stringtie.merged.gtf mergelist.txt &>stringtie.merged.log

cat augustus.gff snap.gff > gene_predictions.gff
perl ./releated_scripts/denovo_change_2_gff3.pl gene_predictions.gff >gene_predictions.gff3

cat genomethreader.gff exonerate.gff > protein_alignments.gff
perl ./releated_scripts/homolog_change_2_gff3.pl protein_alignments.gff >protein_alignments.gff3

perl gff3_gene_prediction_file_validator.pl your.gff3

perl ./EVidenceModeler-1.1.1/EvmUtils/partition_EVM_inputs.pl --genome genome.fasta --gene_predictions gff/gene_predictions.gff3 --protein_alignments gff/protein_alignments.gff3 --transcript_alignments gff/transcripts.fasta.transdecoder.genome.gff3 --segmentSize 100000 --overlapSize 10000 --partition_listing partitions_list.out

perl ./EVidenceModeler-1.1.1/EvmUtils/write_EVM_commands.pl --genome genome.fasta --weights /gpfs/home/yexinhai/project/16_White_Pavo_genome/02_annotation/08_EVM/weights.txt --gene_predictions gff/gene_predictions.gff3 --protein_alignments gff/protein_alignments.gff3 --transcript_alignments gff/transcripts.fasta.transdecoder.genome.gff3 --output_file_name evm.out  --partitions partitions_list.out > commands.list

perl ./EVidenceModeler-1.1.1/EvmUtils/execute_EVM_commands.pl commands.list | tee run.log

bash commands.list

perl ./EVidenceModeler-1.1.1/EvmUtils/recombine_EVM_partial_outputs.pl --partitions partitions_list.out --output_file_name evm.out


perl ./EVidenceModeler-1.1.1/EvmUtils/convert_EVM_outputs_to_GFF3.pl  --partitions partitions_list.out --output evm.out  --genome genome.fasta

find ./ -mindepth 2 -maxdepth 2 -type f -name "evm.out.gff3" -exec cat {} \; >evm.gff3

python ./releated_scripts/gffrename.py evm.gff3 White_Pavo > White_Pavo.gff3

perl ./releated_scripts/getGene.pl --posformat gff White_Pavo.gff3 White_Pavo.genome.fasta > White_Pavo.cds.fasta

perl ./releated_scripts/cds2aa.pl White_Pavo.cds.fasta > White_Pavo.pep.fasta

blastp -query White_Pavo.pep.fasta -db ./uniprot_20180118/uniprot_sprot.20180118.fasta -num_threads 28 -evalue 1e-3 -outfmt 6 -max_target_seqs 1 -out White_Pavo.fasta.blastpout6 &>blastp_run.log

perl ./releated_scripts/add_swissprot_annotation.pl ./uniprot_20180118/uniprot_sprot.20180118.fasta White_Pavo.pep.fasta White_Pavo.fasta.blastpout6  >White_Pavo.pep.anno.fasta
