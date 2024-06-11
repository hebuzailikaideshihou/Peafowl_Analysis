#!/bin/sh
export ref=./White_Pavo.fa
export fastq=./01.all.list
export BAM=./mapping/02.bam
export GVCF=./mapping/03.gvcf

/usr/bin/java -Xmx25g -Djava.io.tmpdir=./tmp -jar $trimmomatic PE -threads 12 \
-summary $fastq/${sample}.summary \
$fastq/${sample}_1.fastq.gz \
$fastq/${sample}_2.fastq.gz \
$fastq/${sample}_1_trimmed.fq.gz \
$fastq/${sample}_1_singleton.fq.gz \
$fastq/${sample}_2_trimmed.fq.gz \
$fastq/${sample}_2_singleton.fq.gz \
LEADING:20 TRAILING:20 SLIDINGWINDOW:3:15 AVGQUAL:20 MINLEN:35 TOPHRED33
echo "trim is Done!"

rm -rf $fastq/${sample}_1_singleton.fq.gz $fastq/${sample}_2_singleton.fq.gz 

bwa mem -t 12 -M -R '@RG\tID:${sample}\tLB:${sample}\tPL:ILLUMINA\tSM:${sample}' $ref $fastq/${sample}_1_trimmed.fq.gz $fastq/${sample}_2_trimmed.fq.gz | ./software/samtools-1.16.1/bin/samtools view -bS - >$BAM/${sample}.bam
echo "mapping is Done!"
/usr/bin/java -Djava.io.tmpdir=./tmp -Xmx50g -jar $PICARD SortSam I=$BAM/${sample}.bam O=$BAM/${sample}.sort.bam SORT_ORDER=coordinate
echo "sort is Done!"
/usr/bin/java -Djava.io.tmpdir=./tmp -Xmx50g -jar $PICARD MarkDuplicates I=$BAM/${sample}.sort.bam O=$BAM/${sample}.sort.dedup.bam M=$BAM/${sample}.marked_dup_metrics.txt REMOVE_DUPLICATES=true CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT
echo "DeDup is Done!"

/usr/bin/java -Djava.io.tmpdir=./tmp -Xmx50g -jar $GATK3_8 \
           -T RealignerTargetCreator \
		   -R $ref \
		   -I $BAM/${sample}.sort.dedup.bam \
		   -nt 12 \
		   -o $BAM/${sample}.RTC.intervals \
		   -allowPotentiallyMisencodedQuals
/usr/bin/java -Djava.io.tmpdir=./tmp -Xmx50g -jar $GATK3_8 \
           -T IndelRealigner \
		   -R $ref \
		   -I $BAM/${sample}.sort.dedup.bam \
		   --targetIntervals $BAM/${sample}.RTC.intervals \
		   -allowPotentiallyMisencodedQuals \
		   -o $BAM/${sample}.realign.bam 

rm -rf $BAM/${sample}.bam $BAM/${sample}.sort.dedup.bam $BAM/${sample}.marked_dup_metrics.txt $BAM/${sample}.RTC.intervals $BAM/${sample}.sort.dedup.bai
            
/usr/bin/java -Djava.io.tmpdir=./tmp -Xmx50g -jar $GATK3_8 \
		-T HaplotypeCaller \
		-R $ref \
		-ERC GVCF \
		-variant_index_type LINEAR \
		-variant_index_parameter 128000 \
		-nct 12 \
		-I $BAM/${sample}.realign.bam \
		-o $GVCF/${sample}.g.vcf.gz

export ref=./White_Pavo.fa
export GVCF=./mapping/03.gvcf
export Genotype=./mapping/04.raw_vcf/79.merge_vcf


/usr/bin/java -Djava.io.tmpdir=./tmp -Xmx100g -jar $GATK3_8 \
        -R $ref \
        -T GenotypeGVCFs \
        --variant $Genotype/78.list \
	-o $Genotype/78.vcf.gz \
	--disable_auto_index_creation_and_locking_when_reading_rods \
	--includeNonVariantSites \
	-nt 12 \
	-stand_call_conf 10

###bamqc

qualimap bamqc -bam $BAM/${sample}.realign.bam -outdir $BAMQC/${sample} -outformat HTML -nt 12 --java-mem-size=50G

###phase

for i in {1..36} MT Z W 
do
	echo $i
	echo "#!/bin/bash" > $i.phase.sh
	echo "export beagle=./software/beagle.27Jan18.7e1.jar" >> $i.phase.sh
	echo "export Finalsnp=./76.snp" >> $i.phase.sh
	echo "java -Xmx40g -Xss128m -jar \$beagle chrom=$i gtgl=\$Finalsnp/$i.clean.SNP.vcf.gz out=chr$i.imp gprobs=true niterations=10 nthreads=6" >> $i.phase.sh
	echo "java -Xmx40g -Xss128m -jar \$beagle gt=chr${i}.imp.vcf.gz out=chr${i}.imp.phase gprobs=true niterations=10 nthreads=6 ibdtrim=40 ibd=true ibdlod=1" >> $i.phase.sh
done

##snpeff

zcat ./all.77.clean.SNP.vcf.gz | grep -v '##' |cut -f 1-9 > White_Pavo.77.header.loc

java -Xmx8g -jar ./software/snpEff/snpEff.jar Pavo White_Pavo.77.header.loc > White_Pavo.77.header.eff.vcf -csvStats White_Pavo.77.header.eff.vcf.csv -stats White_Pavo.77.header.eff.vcf.html

bgzip -c White_Pavo.77.header.eff.vcf > White_Pavo.77.header.eff.vcf.gz

zcat White_Pavo.77.header.eff.vcf.gz | perl vcfEffOnePerLine.pl | java -jar SnpSift.jar extractFields - CHROM POS REF ALT ANN[*].EFFECT ANN[*].IMPACT ANN[*].GENE ANN[*].GENEID ANN[*].HGVS_C ANN[*].HGVS_P ANN[*].RANK ANN[*].DISTANCE > White_Pavo.77.header.eff.split.vcf

bgzip -c White_Pavo.77.header.eff.split.vcf > White_Pavo.77.header.eff.split.vcf.gz

zcat White_Pavo.77.header.eff.split.vcf.gz |grep -E 'MODERATE|HIGH' | bgzip -c > White_Pavo.77.header.eff.split.M_H.vcf.gz
