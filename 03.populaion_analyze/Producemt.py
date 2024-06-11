"""
Created on Sun May 15 20:29:23 2016

"""
def load_name_list(name_file):
	nameList = []
	with open(name_file) as f:
		for name in f:
			nameList.append(name.strip())
		return nameList


def produce(nameList):
	name_num = len(nameList)
	for i in range(name_num):
		with open(str(i+1)+'.SMD.sh','w') as f:
			f.write('''#!/bin/sh
export BAM=
export MTbam=
export mtfq=
export fqmia=
''')
			f.write('''#step1 get mt
samtools view -b ${BAM}/%s.realign.bam MT -h > ${MTbam}/%s.mtgenome.bam
'''%(nameList[i],nameList[i]))
			f.write('''#step2 bam to fq
samtools view -q 10 -b -F 4 ${MTbam}/%s.mtgenome.bam | /usr/bin/java -jar $PICARD SamToFastq INPUT=/dev/stdin FASTQ=${mtfq}/%s.fq INCLUDE_NON_PRIMARY_ALIGNMENTS=false INTERLEAVE=true VALIDATION_STRINGENCY=LENIENT
'''%(nameList[i],nameList[i]))
			f.write('''#step3 assemble MT with MIA
~/software/bin/mia -H 1 -F -i -c -r White_Pavo_MT.fa -f ${mtfq}/%s.fq -m ${fqmia}/%s
'''%(nameList[i],nameList[i]))
			f.close()


def main():
	import sys
	nameList = load_name_list(sys.argv[1])
	produce(nameList)
main()
