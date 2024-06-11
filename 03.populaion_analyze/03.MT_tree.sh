#!/bin/bash
mkdir -p 01_get_mtbam 02.bam_to_fq 03.mia 
cd 01_get_mtbam 
python Producemt.py 78.list
bash *.SMD.sh
cd ..
#step4 to fa 
mkdir -p 04.produce_fa
ls 03.mia > ./04.produce_fa/all.list
cd ./04.produce_fa
ln -s ../03.mia/* .
for i in `cat all.list`; do echo "$i" > list_$i;done 
for i in `cat all.list` ; do j=${i%%.*}; perl ../mtgenome_together.pl list_${i} ${j}.txt; done
rm -rf list*
cat *.txt > 77.raw.fasta
awk -F'.' '{print $1}' 77.raw.fasta > 78.fa
#step5 muscle
mkdir 05.muscle
cd 05.muscle
ln -s ../04.produce_fa/77.fa .
~/software/bin/muscle -in 77.fa -out 78.aln.fa
#step6 tree

seqkit replace -s -p "N" -r "-" 77.aln.fa > 78.alnr.fa

seqkit seq 77.alnr.fa -w 0 > 77.alnr0.fa

~/software/bin/trimal -in 77.alnr0.fa -out 78.trimML.fa -fasta -nogaps
#!/bin/bash
mkdir -p ../05.iqtree

~/software/bin/iqtree2 -s 77.trimML.fa -m MF -nt 8

~/software/bin/iqtree2 -s 77.trimML.fa -m TIM3+F+R2 -nt 8 -bb 1000 

python2 ../seq.file.converter.py -i 77.trimML.fa -inf FASTA -outf PHYLIP > 78.phy

java -jar /public/home/04021/software/jmodeltest-2.1.10/jModelTest.jar -d 77.phy -i -f -g 4 -BIC -AIC -AICc -DT -v -a -w

raxmlHPC-PTHREADS-AVX2 -f a -x 123 -p 23 -# 100 -k -s 77.phy -m GTRGAMMA -n Mit_tree_boot -T 20 -o ERR5432909,ERR5432257