import argparse
import random
from random  import choice
import matplotlib.pyplot as plt
import collections
import pandas as pd


class Population:
    def __init__(self,ne):
        self.ne=ne
        
        if self.ne<=0:
            raise ValueError("Number of diploid individuals should be larger than 0")

    
    @classmethod
    def make_base_pop(cls,pop_size,s_counts): 
        basehaps=[0 for i in range(0,pop_size-s_counts)]+[1 for j in range(0,s_counts)]
        
        basepop=random.sample(basehaps,len(basehaps))
        
        diploids=make_diploids(basepop).convert_to_diploids()
        return(diploids)
        

        
class make_diploids:
    def __init__(self,haplotypes):
        self.haplotypes=haplotypes
        
    def convert_to_diploids(self):    
        x=self.haplotypes[::2]
        y=self.haplotypes[1::2]
        diploid_basepop=zip(x,y)
        
        return(diploid_basepop)

    


class Selection():
    def __init__(self,diploids,s,e):
        
        self.diploids=list(diploids)
        self.s=s
        self.e=e
        self.max_fitness=max(1.0, 1+self.e*self.s, 1+self.s)

    def make_selected_pop(self):
        selected_pop=[]
        
        for i in range(0,len(self.diploids)):
            #print self.diploids[i]
            
            decision=random.random()
            #print decision
            fitness=sum(i for i in self.diploids[i])
            if fitness==2:
                if decision< (1+self.s)/self.max_fitness:
                    selected_pop.append(self.diploids[i])
                else:
                    pass
            elif fitness==1:
                
                if decision < (1+self.s*self.e)/self.max_fitness:
                    selected_pop.append(self.diploids[i])
                else:
                    pass
            elif fitness==0:
                #print fitness
                if decision< 1/self.max_fitness:
                    selected_pop.append(self.diploids[i])
                else:
                    pass
            else:
                raise ValueError("wrong estimation of fitness")

        return(selected_pop)



class mating:
    def __init__(self,selected_pop,ne):
        
        self.selected_pop=list(selected_pop)
        self.ne=ne
        
    def make_mating_pairs(self):
        evolved_pop=[]
 
        while len(evolved_pop) < self.ne:
            ind1,ind2=choice(self.selected_pop),choice(self.selected_pop)
            
           # print ind1,ind2
            hap1=choice(ind1)
            hap2=choice(ind2)
      
            evolved_pop.append((hap1,hap2))
            
        self.selected_pop=evolved_pop
        return(evolved_pop)
    
    def get_frequency(self):
        total_sum=0
        
        for j in range(0, len(self.selected_pop)):
           total_sum+=sum(i for i in self.selected_pop[j])
        
        freq=total_sum/float(len(self.selected_pop*2))
        return(freq)
    



def plotlines(frequencies):
    plt.figure(figsize=(10,7))
    plt.subplots_adjust(left=0.1,bottom=0.1, right=0.9, top=0.9, wspace=0.3, hspace=0.4)

    dict=collections.defaultdict(lambda:[])
    for i in range(0, args.gen+1):
       dict['gen'+'{0}'.format(i)]=frequencies[i::args.gen+1]
    df=pd.DataFrame(dict)
    df.rename(index=lambda x: "repl"+str(x+1), inplace=True)


    mean=df.mean(axis=0)
    variance=df.var(axis=0)
    dfh=(df*2*(1-df)).mean(axis=0)    #df.applymap(lambda x: 2*x*(1-x))

    df.loc["mean",df.columns]=mean
    df.loc["variance",df.columns]=variance
    df.loc["heterozygosity",df.columns]=dfh

    df=df.T
    
    plt.subplot(2,2,1)
    plt.ylim(0,1.0)
    plt.title("Trajectories of random alleles")        
    for i in range(0,len(frequencies), args.gen+1):
        plt.plot(range(0,args.gen+1), frequencies[i:args.gen+1+i],linewidth=0.8 )
        plt.ylabel("Frequency")

    plt.subplot(2,2,2)
    plt.title("Mean allele frequency")
    plt.plot(range(0,args.gen+1),df["mean"],color="black",linestyle="--", linewidth=2.0)
    plt.ylim(0,1.0)

    plt.subplot(2,2,3)
    plt.title("Variance among replicates")
    plt.plot(range(0,args.gen+1),df["variance"],color="black",linestyle="--", linewidth=2.0)
    plt.ylim(0,0.5)
    plt.xlabel("Generations")
    plt.ylabel("Variance")


    plt.subplot(2,2,4)
    plt.title("Mean Heterozygosity")
    plt.plot(range(0,args.gen+1),df["heterozygosity"],color="black",linestyle="--", linewidth=2.0)
    plt.ylim(0,0.6)
    plt.xlabel("Generations")
    plt.ylabel("Heterozygosity")


    
    plt.savefig("../results/{0}Ne_{1}freq_{2}gen_{3}repl.png".format(args.ne,args.p,args.gen,args.repl))


 

       
def plotresults(freqs):
    bins=75
    ydata=frequencies[args.gen::args.gen+1]
    print(ydata)
    ydata= map(lambda x: x*args.ne*2, ydata)
    ydata=list(map(int,ydata))
    print(ydata)
    plt.hist(ydata, bins, histtype='bar', rwidth=4)
    plt.show()

    

parser=argparse.ArgumentParser(description= """
            Description
            -----------
            Python script that simulates the evolutionary trajectory of a locus (for diploids).
            
            Authors 
            -----------
            Vlachos Christos""",formatter_class=argparse.RawDescriptionHelpFormatter)


parser.add_argument("-Ne",type=int, required=True,dest="ne",default=None, help="Population Size")
parser.add_argument("-s",type=float, required=True,dest="s",default=None, help="Selection Coefficient")
parser.add_argument("-p",type=float, required=True,dest="p",default=None, help="Initial Allele Frequency")
parser.add_argument("-e",type=float, required=True,dest="e",default=None, help="Dominance effect")
parser.add_argument("--replicates",type=int, required=True,dest="repl",default=None, help="Number of replicates")
parser.add_argument("--generations",type=int, required=True,dest="gen",default=None, help="Number of generations")
parser.add_argument("--drift-distribution", required=False,dest="distr",default=False, help="If True returns the distribution of allele frequencies in the last generation (s must be 0)")

args = parser.parse_args()

pop_size=2*args.ne
s_counts=int(args.p*pop_size)
frequencies=[]


for i in range(0,args.repl):
     population=list(Population(args.ne).make_base_pop(pop_size, s_counts))
#     print(population)
     freq=mating(population,args.ne).get_frequency()
     frequencies.append(freq)
    
     for j in range(0,args.gen):
         print("Processing generation {0} of replicate {1}".format(j+1,i+1))
         population=Selection(population,args.s,args.e).make_selected_pop()
#         print(population)
         population=mating(population,args.ne).make_mating_pairs()
#         print(population)
         freq=mating(population,args.ne).get_frequency()

         frequencies.append(freq)
 
#print len(frequencies)
#print frequencies
          
if args.s==0.0 and args.distr=="True": 
    plotresults(frequencies)    
else:  
    plotlines(frequencies)
    

    
    
    
    
    
    
    
