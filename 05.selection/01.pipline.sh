#!/bin/bash
vcftools --gzvcf ./all.18.clean.new.SNP.vcf.gz \
	--weir-fst-pop ./BP.list \
	--weir-fst-pop ./WP.list \
	--fst-window-size 1000 \
	--fst-window-step 200 \
	--out Fst-win \
	--max-missing 0.9 \
	--maf 0.05

bcftools filter ./4.clean.SNP.vcf.gz  --regions 4 > EDNRB2.vcf

awk 'BEGIN { FS = "\t"; OFS = "\t" } { for (i=1; i<=NF; i++) { if ($i ~ /^\.:.*/) $i = "./.:" substr($i, 3) } }1' EDNRB2.vcf >EDNRB2_corr.vcf

python parseVCF.py -i  EDNRB2_corr.vcf | bgzip > introgression.geno.gz


python popgenWindows.py -w 1000 -s 200 -g introgression.geno.gz -o output.csv.gz -f phased -T 5 -p popA -p popB -p popC --popsFile pops.txt

###RNA_seq
ref=./White_Pavo.fa
gff=./genes.gtf
index=./White_Pavo
label=
out=
mkdir -p ${out}
fq=
fastp -i ${fq}/${label}.fq.gz -o ${fq}/${label}.clean.fq.gz
hisat2 --dta -p 6 -U ${fq}/${label}.clean.fq.gz -x ${index} -S ${out}.sam
samtools view -@ 6 -b -S ${out}.sam > ${out}.bam
samtools sort -@ 6 -m 4G -o ${out}.sorted.bam ${out}.bam
stringtie -p 6 -G ${gff}  -l ${label} -o ${out}/${out}.gtf ${out}.sorted.bam
stringtie -e -B -p 6 ${i}.sorted.bam -G stringtie_merged.gtf -o ./Ballgown/${i}/${i}.gtf
ls ./Ballgown/*/*.gtf > ./8.gtf.list
awk -F "/" '{print $3"\t"$0}' ./8.gtf.list >./8.sample.gtf.list
python2 ./prepDE.py -i 8.sample.gtf.list
Rscript DEseq2.R
