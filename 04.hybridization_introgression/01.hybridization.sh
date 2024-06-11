#!/bin/bash
###Identification of purebred peacocks using the results of f3, Dsuite and admixture
qp3Pop -p pop3.3.par > f3.3.out
###loter

for i in {1..36} Z
do
    python ./01_loter_multiRegions.py  ./77.imp.phase.vcf.gz  ${i}.hy --groupfile 77.pop --regionfile ${i}.region --refpops WP,GP --querypops GFf1A,qhd01,qhd02 --threads 8 --nbags 20
done

for i in {1..36} Z
do
    python ./02_merge_loter_multiRegions.py ${i}.hy.anc.tsv.gz ${i}.indseg.tsv.gz
done

for i in {1..36} Z
do
    python ./03_add_meanBag_to_indSeg.py ${i}.indseg.tsv.gz ${i}.hy.bag.tsv.bgz ${i}.indseg.addBag.tsv.gz
done

for i in {1..36} Z
do
    python ./04_cal_freq_persite.py --ancfile ${i}.hy.anc.tsv.gz --bagfile ${i}.hy.bag.tsv.gz --outfile ${i}.hy.freq.tsv.gz --popfile hy.pop --chunksize 10000 --cutoff 144 --refpops WP,GP
done  
python ./04_cal_freq_persite.py --ancfile 1.hy.anc.tsv.gz --bagfile 1.hy.bag.tsv.gz --outfile 1.hy.freq.tsv.gz --popfile hy.pop --chunksize 10000 --cutoff 144 --refpops WP,GP

###Dsuite
module load compiler/gcc/7.3.1
./Dsuite-master/Build/Dsuite Dtriosall.78.clean.SNP.vcf.gz 78.pop
./Dsuite-master/Build/Dsuite Dinvestigate -w 2500,500 ./all.78.clean.SNP.vcf.gz 78.pop test_trios.txt

###plot
./RectChr-1.34/Example/example4/../../bin/RectChr    -InConfi    in.cofi   -OutPut OUT