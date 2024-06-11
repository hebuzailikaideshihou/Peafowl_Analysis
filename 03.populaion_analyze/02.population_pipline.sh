#!/bin/bash


####  Vcf to  plink format
vcftools --gzvcf ./all.77.clean.SNP.vcf.gz \
	--keep kongque-56.list \
        --plink \
	--chrom-map ChrID \
	--out kongque-56
  
#### make bed
  
plink --file ./kongque-56 --make-bed --chr-set 36 --out kongque-56 
  
#### Remove LD sites

plink -bfile ./kongque-56 --indep-pairwise 50 5 0.2 --chr-set 36 --out kongque-56
   
## extract non_LD sites
plink -bfile ./kongque-56 --extract kongque-56.prune.in --make-bed --chr-set 36 --out kongque-56.prune
  
### Construct NJtree

plink -bfile ./kongque-56 --chr-set 36 --distance-matrix --out kongque-56
 
perl plink.distance.matrix.to.mega.pl kongque-56.mdist.id kongque-56.mdist 56 kongque-56.NJtree

for i in {1..10}
do
	admixture --cv ./kongque-56.prune.bed $i | tee log${i}.out
   	admixture ./kongque-56.prune.bed $i
done

grep -h CV log*.out > CV.value
paste kongque-56.prune.2.Q kongque-56.prune.3.Q kongque-56.prune.4.Q kongque-56.prune.5.Q > kongque-56.prune.Q
## PCA
smartpca -p Kongque_eigensoft_smartpca.par
less kongque-56.SNP.evec  |awk -F ':' '{print $1}'| sed 's/\ //g' > sample.txt
paste sample.txt kongque-56.SNP.evec | grep -v "#"  > kongque-56.SNP.pca 
VCF2PCACluster \
        -InVCF ./all.77.clean.SNP.vcf.gz \
        -OutPut kongque-56  \
        -MAF 0.001 \
        -Fchr Z,W,MT \
        -InSampleGroup  kongque-56.pop

###ROH
#!/bin/bash 
for i in BPB BPJ BPW GQHD GTL GYN GYUEN GZJ WPB WPW 
do
	vcftools --gzvcf /public/home/04021/01.chinese/01.wg/01.kongque/mapping/05.filter/77.snp/all.77.clean.SNP.vcf.gz --not-chr MT --not-chr Z --not-chr W --keep ${i}.list --plink --chrom-map ChrID --out ${i}_snp
	plink --file ${i}_snp --make-bed --chr-set 36 --out ${i}_snp
	plink --bfile ${i}_snp --chr-set 36 --homozyg-window-snp 50 --homozyg-snp 100 --homozyg-kb 500 --homozyg-density 50 --homozyg-gap 100 --homozyg-window-missing 2 --homozyg-window-threshold 0.05 --homozyg-window-het 1 -out ${i}_snp 
done

#!/bin/bash 
#for i in  `cat ./breed `
#do 
#	grep "$i" ~/Donkey/01.selection-signal/02.liangzhou/00.pop/1.Liangzhou-59.list | awk '{print $3}' > $i.list
#done

for i in `cat ./breed `
do 
	cat ${i}_snp.hom | grep -v 'FID' | awk '{print $1,$9}' >> all.hom.1
done

sed -i '1i FID\tKB' all.hom.1

for i in `cat ./breed `
do
	cat ${i}_snp.hom | grep -v 'FID' > ${i}_snp.hom1
	python ./count_length_roh.1.py ${i}_snp.hom1 ${i}.sum
	rm -rf ${i}_snp.hom1
done
for i in `cat ./breed `
do
	cat ${i}_snp.hom.indiv | grep -v 'FID' >> all.indiv
done

cat WP_snp.hom.indiv | grep 'FID' > indiv.head

cat indiv.head all.indiv > all.indiv.1
#!/bin/bash

for i in `cat 68.list`
do cat all.hom.1 |awk '$1=="'$i'" {print $0}' >> $i.hom
done

for i in `cat 68.list` 
do cat $i.hom |awk 'NR==1{min=$2;next}{min=min<$2?min:$2}END{print $i "\t" min}' > $i.min
done

for i in `cat 68.list` 
do cat $i.hom |awk 'BEGIN {max = 0} {if ($2+0 > max+0) max=$2} END {print $i,"Max=", max}' > $i.max
done

for i in `cat 68.list`
do paste $i.min >> all.min
done

for i in `cat 68.list`
do paste $i.max >> all.max
done

for i in `cat 68.list`
do
	cat $i.hom | awk '{sum+=$2} END {print "'${i}'",sum,sum/NR, sum/2432132968}' >>  all.roh_sum
done

###LD
###STEP1
for i in `cat pop.list`
do 
    PopLDdecay -InVCF ./all.77.clean.SNP.vcf.gz -SubPop ${i}.list  -OutStat ${i}.res
done
#!/bin/bash
###STEP2
for i in `cat pop.list`
do
    echo "`pwd`/$i.res.stat.gz  ${i}" >> plot.list
done
echo "Plot_MultiPop.pl -inList plot.list  -output all" >> plot.sh
chmod 755 plot.sh
bash plot.sh

###Pi
#!/bin/bash 
for i in `cat pop.list`
do
	vcftools --gzvcf ./all.77.clean.SNP.vcf.gz \
		--keep ${i}.list
        --window-pi 50000 \
		--out ${i}
	grep -v "CHR" ${i}.windowed.pi | awk '{print "'"${i}"'""\t"$1"\t"$2"\t"$3"\t"$4"\t"$5}'> ${i}.pi
	less ${i}.windowed.pi | grep -v "CHR" | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"-(log($5)/log(10))}' > ${i}.ratio.pi
	awk '{print "'"${i}"'""\t"$1"\t"$2"\t"$3"\t"$4"\t"$5}' ${i}.ratio.pi > ${i}.ratio_breed.pi
done

cat header *.ratio.pi > all.ratio_breed.pi
cat header *.pi > all.pi