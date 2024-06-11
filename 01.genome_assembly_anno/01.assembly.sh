#!/bin/bash
## assembly software: hifiasm
## Hic-anchor software: yahs

################################
#- workdir
#	- Hifi_input/
#	- HiC_input/	
#	- Hifi_output/
#	- juicer_output/ 
################################

##shell
while getopts ":f:c:t:" opt
do
	case $opt in
		f)
			Hifi_input=$OPTARG;;
		c)
			HiC_input=$OPTARG;;
		t)
			CPU=$OPTARG;;
		?)
			echo "Usage: assembly_anchor.sh -f hifi_input -c hic_input -t cpu"
			exit 1;;
	esac
done

#Hifi_input=$1
#HiC_input=$2
#CPU=$3

workdir="/public/home/04021/01.chinese/01.wg/01.kongque/Denove/06.test"
Prefix="White_Pavo"

juicer="/public/home/04021/software/juicer-1.6"
ddna3="/public/home/04021/software/3d-dna-201008"

# Assembly
## create Assembly_OUTPUT directory
if [ ! -d "${workdir}/Assembly_OUTPUT" ];then
	mkdir ${workdir}/Assembly_OUTPUT
else
	break
fi

hifiasm -o ${workdir}/Assembly_OUTPUT/$Prefix -t $CPU $(echo ${workdir}/${Hifi_input}/*) && echo "**** hifiasm done! ****"
awk '/^S/{print">"$2;print $3}' ${workdir}/Assembly_OUTPUT/$Prefix.bp.p_ctg.gfa > ${workdir}/Assembly_OUTPUT/$Prefix.bp.p_ctg.fasta && echo "**** ${Prefix}.bp.p_ctg.fasta done! ***"

#Index 
bwa index ${workdir}/Assembly_OUTPUT/$Prefix.bp.p_ctg.fasta && echo "**** bwa index done! ****"

#Hic-anchor
# juicer

## create juicer_output & data directory and move Hi-C data to data directory
if [ ! -d "${workdir}/juicer_output" ];then
	mkdir -p ${workdir}/juicer_output/data/fastq
	ln -s ${workdir}/${HiC_input}/* ${workdir}/juicer_output/data/fastq
else
	break
fi

### contig
awk 'BEGIN{OFS="\t"}{print $1, $NF}' ${workdir}/juicer_output/${Prefix}_HindIII.txt > ${workdir}/juicer_output/${Prefix}.chrom.sizes && echo "**** ${Prefix}.chrom.sizes done! ****"

## need juicer_tools/pretextmap and samtools if want to do hic plot
## juicer_tools: https://github.com/aidenlab/juicer/wiki/Download
## PretextMap: https://github.com/wtsi-hpag/PretextMap
## PretexSnapshot: https://github.com/wtsi-hpag/PretextSnapshot
## samtools: https://github.com/samtools/samtools
## please adjust the path to juicer_tools and samtools
## here we use 12 CPUs and 32Gb memory for juicer_tools pre - adjust it according to your device
## see more information for juicer tools https://github.com/aidenlab/juicer/wiki/Juicer-Tools-Quick-Start
## output file will be ${outdir}/${out}.hic
## the output hic file could be viewed with JuiceBox https://github.com/aidenlab/Juicebox
noplot=0
#juicer_tools="java -Xmx32G -jar /bin/juicer_tools_1.22.01.jar pre --threads 12"
## v1.9.9 seems much faster than v1.22.01
juicer_tools="java -Xmx400G -jar /public/home/04021/01.chinese/01.wg/software/juicer_tools_1.19.02.jar pre"
pretext_map="/public/home/04021/anaconda3/envs/whatshap/bin/PretextMap"
pretext_snapshot="/public/home/04021/anaconda3/envs/whatshap/bin/PretextSnapshot"
samtools="/public/home/04021/software/samtools-1.16.1/bin/samtools"

#### run yahs scaffolding
#yahs -r 1000,2000,5000,10000,20000,50000,100000,200000,500000 -o ${outdir}/${out} ${contigs} ${hicaln} >${outdir}/${out}.log 2>&1
yahs -e AAGCTAGCTT -o ${outdir}/${out} ${contigs} ${hicaln} >${outdir}/${out}.log 2>&1

if [ ${noplot} -ne 0 ]; then exit 0; fi

#### this is to generate input file for juicer_tools - non-assembly mode or for PretextMap
## here we use 8 CPUs and 32Gb memory for sorting - adjust it according to your device
(juicer pre ${outdir}/${out}.bin ${outdir}/${out}_scaffolds_final.agp ${contigs}.fai 2>${outdir}/tmp_juicer_pre.log | LC_ALL=C sort -k2,2d -k6,6d -T ${outdir} --parallel=20 -S400G | awk 'NF' > ${outdir}/alignments_sorted.txt.part) && (mv ${outdir}/alignments_sorted.txt.part ${outdir}/alignments_sorted.txt)
## prepare chromosome size file from samtools index file
# ${samtools} faidx ${outdir}/${out}_scaffolds_final.fa
# cut -f1-2 ${outdir}/${out}_scaffolds_final.fa.fai >${outdir}/${out}_scaffolds_final.chrom.sizes
## another way to prepare chromosome size file
## this is an easier way especially when we have >2G scaffolds which need scaling 
cat ${outdir}/tmp_juicer_pre.log | grep "PRE_C_SIZE" | cut -d' ' -f2- >${outdir}/${out}_scaffolds_final.chrom.sizes
## do juicer hic map
(${juicer_tools} ${outdir}/alignments_sorted.txt ${outdir}/${out}.hic.part ${outdir}/${out}_scaffolds_final.chrom.sizes) && (mv ${outdir}/${out}.hic.part ${outdir}/${out}.hic)
## do Pretext hic map
(awk 'BEGIN{print "## pairs format v1.0"} {print "#chromsize:\t"$1"\t"$2} END {print "#columns:\treadID\tchr1\tpos1\tchr2\tpos2\tstrand1\tstrand2"}' ${outdir}/${out}_scaffolds_final.chrom.sizes; awk '{print ".\t"$2"\t"$3"\t"$6"\t"$7"\t.\t."}' alignments_sorted.txt) | ${pretext_map} -o ${outdir}/${out}.pretext
# and a pretext snapshot
${pretext_snapshot} -m ${outdir}/${out}.pretext --sequences "=full" -o ${outdir}

#### this is to generate input file for juicer_tools - assembly (JBAT) mode (-a)
juicer pre -a -o ${outdir}/${out}_JBAT ${outdir}/${out}.bin ${outdir}/${out}_scaffolds_final.agp ${contigs}.fai 2>${outdir}/tmp_juicer_pre_JBAT.log
cat ${outdir}/tmp_juicer_pre_JBAT.log | grep "PRE_C_SIZE" | cut -d' ' -f2- >${outdir}/${out}_JBAT.chrom.sizes
(${juicer_tools} ${outdir}/${out}_JBAT.txt ${outdir}/${out}_JBAT.hic.part ${outdir}/${out}_JBAT.chrom.sizes) && (mv ${outdir}/${out}_JBAT.hic.part ${outdir}/${out}_JBAT.hic)

#### this is to generate final genome assembly file after manual curation with JuiceBox (JBAT)
## the output assembly file after curation is ${outdir}/${out}_JBAT.review.assembly
## the final output is ${outdir}/${out}_JBAT.FINAL.agp and ${outdir}/${out}_JBAT.FINAL.fa
juicer post -o ${outdir}/${out}_JBAT ${outdir}/${out}_JBAT.review.assembly ${outdir}/${out}_JBAT.liftover.agp ${contigs}