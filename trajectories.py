import argparse
import gzip
import random
import matplotlib.pyplot as plt



def selected_line(sel_line):

    k=sel_line.split("\t")

    chrm_sel=k[0]
    pos_sel=k[1]
    #selected_minus,selected_plus=k[2].split("/")
    selected_plus,selected_minus=k[2].split("/")
    effect_size=k[3]
    dom_effect=k[4]

    return (chrm_sel,pos_sel,selected_plus,effect_size)



def all_freqs(line,selected_plus):
    freq_list=[]
    

    k=line.split("\t")
    chrm=k[0]
    pos=k[1]
    snp=k[2]



    freq_dict={}
    for i in range(3, int(rep_number)+3):
        generations_count=k[i]
        count_A_gen,count_T_gen,count_C_gen,count_G_gen,count_N_gen,count_mut_gen= generations_count.split(":")
        total_count=sum([int(count_A_gen),int(count_T_gen),int(count_C_gen),int(count_G_gen),int(count_N_gen),int(count_mut_gen)])

        freq_dict={'A':count_A_gen,'T':count_T_gen,'C':count_C_gen,'G':count_G_gen}
        generations_freq=float(freq_dict.get(selected_plus))/total_count
        freq_list.append(generations_freq)
        
        generations_counts=[]

    return (freq_list)



def mean_freqs(freqs):
	gen_fqs=[]
	rep_fqs=[]
	for rep in freqs:
		for gen in range(0,args.tp):
			fq=rep[int(gen)::args.tp]
			mean_fq=round(sum(fq)/len(fq),3)
			gen_fqs.append(mean_fq)

		rep_fqs.append(gen_fqs)

		gen_fqs=[]

	return(rep_fqs)
	





fig, axes = plt.subplots(nrows=1,ncols=2,figsize=(12,6))
def plotlines(mfqs, pl):
	if pl==0:
		for snp in mfqs:
			axes[pl].title.set_text('Average frequency of selected SNPs')
			axes[pl].set_ylabel('Frequency')
			axes[pl].set_xlabel('Generations')
			plo=axes[pl].plot(list(range(0,rep_number,10)),snp, 'b--')
			 
		return(plo)

	elif pl==1:
		for snp in mfqs:
			axes[pl].set_ylabel('Frequency')
			axes[pl].set_xlabel('Generations')
			axes[pl].title.set_text('Average frequency of neutral SNPs')
			axes[pl].set_ylim([0, 1])
			plo=axes[pl].plot(list(range(0,rep_number,10)),snp, 'black')

		return(plo)

	else:
		raise ValueError("The plot id is wrong")
		return
		





parser=argparse.ArgumentParser(description= """
            Description
            -----------
            Python script that estimates allele frequencies from sync files and returns the average across replicates

            Authors
            -----------
            Vlachos Christos""",formatter_class=argparse.RawDescriptionHelpFormatter)


parser.add_argument("--sync-file", type=str, required=True, dest="sync", default=None, help="MimicreEE2 haplotypes; either *.gz or not")
parser.add_argument("--chosen-alleles",type=str, dest="chosen_alleles", default=None, help="Input file with chosen alleles")
parser.add_argument("--replicates",type=int, dest="repl", default=10, help="Number of replicates in sync file")
parser.add_argument("--time-points",type=int, dest="tp", default=6, help="Number of time-points in sync file")
args = parser.parse_args()



###open file and count columns####
fh=None
myfile=args.sync
if(myfile.endswith(".gz")):
    fh=gzip.open(myfile,mode='rb')
else:
    fh=open(myfile)        


rep_number=args.repl*args.tp



dictionary={}
sel_freqs=[]
neutral_freqs=[]


for sel_line in open(args.chosen_alleles, "r"):
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
       selected_plus=(dictionary[chrm+pos][2])
       sel_freqs.append(all_freqs(line,selected_plus))


    else:
       selected_plus=random.sample(['A','C'],1)[0]
       neutral_freqs.append(all_freqs(line,selected_plus))
       


sel_mfs=mean_freqs(sel_freqs)
neutral_mfs=mean_freqs(neutral_freqs)

s=plotlines(sel_mfs,0)
n=plotlines(neutral_mfs,1)

fig.subplots_adjust(left=0.1,bottom=0.1, right=0.9, top=0.9, wspace=2.9, hspace=1.7)
fig.tight_layout()
plt.show()

