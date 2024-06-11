#!/bin/bash
orthofinder -f ./ -S diamond  -t 24 -a 200

trimal -in SpeciesTreeAlignment.fa -out SpeciesTreeAlignment_trim.fa -fasta -automated1
python2 seq.file.converter.py -i SpeciesTreeAlignment_trim.fa -inf FASTA -outf PHYLIP > SpeciesTreeAlignment_trim.phy
iqtree -s SpeciesTreeAlignment_trim.fa -bb 1000 -bnni -nt 8 -m MFP
##Timetree
baseml baseml.ctl
mcmctree mcmctree3.ctl
mv out.BV in.BV
mcmctree mcmctree2.ctl

###cafe gene family
cafe5 -i cafe.tab -t tree.txt -p -k 2 -o k3p
cafeplotter -i k3p -o k3p_plot

###PSG and RGS

python ./callCodeml.py ./ tree.nwk 64

python ./callCodeml_RGS.py ./ tree.nwk 64

###catus
cactus --workDir=./02.catus --batchSystem single_machine --gpu 1 --lastzCores 2  --maxCores 2  jobStore_foursample seqfile.txt forsample.hal --binariesMode local

###psmc
#!/bin/bash
for i in BP-1B WP-1B
do 
    bcftools mpileup -Ou -I -f White_Pavo.fa ${i}.realign.bam | bcftools call -c -Ov | vcfutils.pl vcf2fq -d 30 -D 200| gzip > ${i}.fq.gz
    ./psmc-master/utils/fq2psmcfa -q20 ${i}.fq.gz > ${i}.psmcfa
    ./psmc-master/psmc -N30 -t5 -r5 -p "4+30*2+4+6+10" -o ${i}.psmc ${i}.psmcfa
done
python  plot_py/easy_psmc_plot.py --psmc_file_list test.list --mutation_rate 1.33e-9 --generation_time 4 --x_max 0.5e7  --y_min 1e2 --y_max 0.3e13 --x_min 1 --span_min 1.5e3 --span_max 1.9e3 --span_color yellow

###Synteny_collinearity
perl ./NGenomeSyn-1.41/bin/GetTwoGenomeSyn_1000.pl    -InGenomeA  chicken.chr.fa -InGenomeB  White_Pavo.fa   -OutPrefix  chicken_White_Pavo