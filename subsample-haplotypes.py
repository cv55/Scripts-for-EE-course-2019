import argparse
import random

parser=argparse.ArgumentParser(description= """
            Description
            -----------
            Python script that picks N neutral SNPs from the genome

            Authors
            -----------
            Vlachos Christos""",formatter_class=argparse.RawDescriptionHelpFormatter)


parser.add_argument("--haplotypes", type=str, required=True, dest="haps", default=None, help="MimicreEE2 haplotypes; either *.gz or not")
parser.add_argument("--chosen-alleles",type=str, dest="selected", default=None, help="Input file with chosen alleles")
parser.add_argument("--output",type=str, dest="out", default=None, help="Name of the output file with the subsampled haplotypes")
args = parser.parse_args()


####open file####
fh=None
myfile=args.haps
if(myfile.endswith(".gz")):
    fh=gzip.open(myfile,mode='rb')
else:
    fh=open(myfile)


writefile=open(args.out, 'w')
dictionary={}
neutral=[]
count_ancestral=0
count_derived=0



def selected_line(sel_line):

    k=sel_line.split("\t")

    chrm_sel=k[0]
    pos_sel=k[1]
    selected_plus,selected_minus=k[2].split("/")
    effect_size=k[3]
    dom_effect=k[4]

    return (chrm_sel,pos_sel,selected_plus,effect_size)



for sel_line in open(args.selected, "r"):
    sel_line=sel_line.rstrip()

    (chrm_sel,pos_sel,selected_plus,effect_size)=selected_line(sel_line)
    key="{0}{1}".format(chrm_sel,pos_sel)
    dictionary[key]=selected_line(sel_line)


for line in fh:
	if(myfile.endswith(".gz")):
		line=line.decode(encoding='utf-8')
	else:
		line=line
		k=line.split("\t")
		chrm=k[0]
		pos=k[1]


	if chrm+pos in dictionary.keys():
		writefile.writelines(line)

	else:
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


		if freq_ancestral>0.75 or  freq_derived>0.75:
			pass
		else:
			neutral.append(line)

		count_ancestral, count_derived= 0,0

    	

neutral=random.sample(neutral, 200)



for i in neutral:
	writefile.writelines(i)

writefile.close()
