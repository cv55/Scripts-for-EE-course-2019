from scipy.stats import norm 
import argparse
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd



parser=argparse.ArgumentParser(description= """
            Description
            -----------
            Python script that estimates the scaled fitness function for stabilizing selection
            
            Authors
            -----------
            Vlachos Christos""",formatter_class=argparse.RawDescriptionHelpFormatter)


parser.add_argument("--qff",type=str, required=True,dest="qff",default=None, help="Quantitative fitness function in the form: minPheno:maxPheno:mean:sd, eg: --qff=-5:5:0:1")
parser.add_argument("--gpf", type=str,required=True, dest="gpf", default=None, help="A GPF output file from MimicrEE2")
args = parser.parse_args()



fig, axes = plt.subplots(nrows=1,ncols=2,figsize=(12,6))
secondary_ax = []


class getFitness:
        def __init__(self,qff):
                

                self.maxFit=1
                self.minFit=0.0
                self.minPhe=qff[0]
                self.maxPhe=qff[1]
                self.mean=qff[2]
                self.sd=qff[3]				
                self.ticks=np.arange(self.minPhe,self.maxPhe, 0.01)

        def estimation(self):
                fig.tight_layout() 
                scl=norm.pdf(self.mean, loc=self.mean, scale=self.sd)
                fitness=norm.pdf(self.ticks, loc=self.mean, scale=self.sd)
                diff=self.maxFit-self.minFit
                sc_fitness=(fitness*diff/scl) +self.minFit
                color = 'tab:red'
                axes[0].set_xlabel('Phenotype')
                axes[0].set_ylabel('fitness', color=color)
                axes[0].tick_params(axis='y', labelcolor=color)
                ff=axes[0].plot(self.ticks,sc_fitness,'r--')
                
                return(ff)


def basepop(gpf):
	axes2 = axes[0].twinx()
	color = 'tab:blue'
	axes2.set_ylabel('Number of individuals', color=color)
	axes2.tick_params(axis='y', labelcolor=color)
	gpf=pd.read_csv(gpf,sep="\t", names=['Replicate', 'Generation','sex', 'Genotype', 'Phenotype', 'Fitness'])
	bp=gpf.head(500)
	pheno=bp['Phenotype']
	bp=axes2.hist(pheno,alpha=0.25)
	return(bp)


def phenotypes(gpf,mean):

	gpf=pd.read_csv(gpf,sep="\t", names=['Replicate', 'Generation','sex', 'Genotype', 'Phenotype', 'Fitness'])
	mean_df=pd.DataFrame(columns=['Generation', 'Replicate', 'Phenotype'], index=range(0,int(gpf.shape[0]/500)))
	c=0
	pheno=[]
	phenos=[]
	gens=[]
	for gen in range(0,60,10):
		gpf_gen=gpf.loc[gpf['Generation'] == gen]
		for rep in range(1,11,1):
			gpf_gen_rep=gpf_gen.loc[gpf['Replicate'] == rep] 

			mean_gen_rep=gpf_gen_rep.Phenotype.mean()
			
			mean_df.loc[c,'Generation']=gen
			mean_df.loc[c,'Replicate']=rep
			mean_df.loc[c,'Phenotype']=mean_gen_rep
			c+=1

			pheno.append(mean_gen_rep)	
		phenos.append(pheno)
		pheno=[]
		gens.append("gen{0}".format(gen) )			

#	print(phenos)
	axes[1].tick_params(labelrotation=45)
	axes[1].set_ylabel('Phenotype')
	ph=axes[1].boxplot(phenos, labels=gens)
	ph=axes[1].axhline(y=mean, linewidth=1,linestyle='--', color = 'black')
		
	
	return(ph)


qff=list(map(float, list(args.qff.split(':'))))
distr=getFitness(qff)
ff=distr.estimation()
bp=basepop(args.gpf)


phenotypes(args.gpf, qff[2])

fig.subplots_adjust(left=0.1,bottom=0.1, right=0.9, top=0.9, wspace=2.9, hspace=1.7)
fig.tight_layout() 
plt.show()



