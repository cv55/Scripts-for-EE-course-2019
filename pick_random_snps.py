import sys
import random
from optparse import OptionParser, OptionGroup
import gzip
import numpy
import timeit

start = timeit.default_timer()

parser = OptionParser()
parser.add_option("--input",dest="input",help="A file containing the dgrp haplotypes (MimicrEE input)")
parser.add_option("--output",dest="output",help="A file containing the randomly selected SNPs")
parser.add_option("--snp-number",dest="snps",type=int, help="Number of selected SNPs")
(options, args) = parser.parse_args()

rndm=None
gamma_list=[]
line_to_parse_freqval_2L=[]
line_to_parse_freqval_3R=[]
count_ancestral=0
count_derived=0
rndm2=None


writefile=open(options.output, 'w')


##keep one line after every ....

for line in open(options.input, 'r'):
    if line.startswith("2L") or line.startswith("2R") or line.startswith("3L") or line.startswith("3R") :
        rndm2=random.random()


        if rndm2<=float(0.01):
            filter_line=line
        #print filter_line


            k=filter_line.split("\t")

            chrm=k[0]
            pos=k[1]
            ancestral,derived=k[3].split("/")
            alleles=k[4].split(" ")


            for j in range(0, len(alleles)):
                allele_to_count=alleles[j]

                if allele_to_count.rstrip()==str(ancestral*2):
                    count_ancestral=count_ancestral+1

                elif allele_to_count.rstrip()==str(derived*2):
                    count_derived=count_derived+1

            freq_ancestral=float(count_ancestral) / len(alleles)
            freq_derived=float(count_derived) / len(alleles)


            if chrm=="2L" and ( freq_ancestral>0.07 and freq_ancestral<0.3 or freq_derived>0.07 and freq_derived<0.3): 
             #  pass
                 
            #else:
                line_to_parse_freqval_2L.append(line)

            elif chrm=="3R" and ( freq_ancestral>0.07 and freq_ancestral<0.3 or freq_derived>0.07 and freq_derived<0.3):
                line_to_parse_freqval_3R.append(line)

            count_ancestral, count_derived= 0,0

##keep random lines(snps)

if len(line_to_parse_freqval_2L and line_to_parse_freqval_3R)>= int(options.snps/2):
    random_lines_2L=random.sample(line_to_parse_freqval_2L, int(options.snps/2))
    random_lines_3R=random.sample(line_to_parse_freqval_3R, int(options.snps/2))

else:
    raise ValueError("Number of random snps greater than the total number of alleles")





def make_file(random_lines, snps):
	count_ancestral, count_derived= 0,0
	for j in range(0, snps):

		line_to_parse2=random_lines[j]
		k2=line_to_parse2.split("\t")

		chrm2=k2[0]
		pos2=k2[1]
		ancestral2,derived2=k2[3].split("/")
		alleles2=k2[4].split(" ")

		for j in range(0, len(alleles2)):
			allele_to_count=alleles2[j]
			if allele_to_count.rstrip()==str(ancestral*2):
				count_ancestral=count_ancestral+1
			elif allele_to_count.rstrip()==str(derived*2):
				count_derived=count_derived+1


		freq_ancestral=float(count_ancestral) / len(alleles)
		freq_derived=float(count_derived) / len(alleles)
		print(freq_ancestral,freq_derived)

		if freq_ancestral>=0.5:
			major=ancestral2
			minor=derived2
		else:
			major=derived2
			minor=ancestral
		effect=1

		writefile.writelines(chrm2+"\t"+ pos2 +"\t"+minor+"/"+major+"\t"+"%f" %effect+"\t"+"%d" %0+"\n")

		count_ancestral, count_derived= 0,0

make_file(random_lines_2L, 5)
make_file(random_lines_3R, 5)


writefile.close()



stop = timeit.default_timer()

print("Elapsed time (sec):", (stop - start))
