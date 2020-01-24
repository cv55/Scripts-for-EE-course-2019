import rpy2.robjects as ro
import argparse
#import pandas as pd
import rpy2.robjects.lib.ggplot2 as ggplot2
from rpy2.robjects.packages import importr
import scipy.special as ss
import math
import numpy as np
import matplotlib.pyplot as plt




parser = argparse.ArgumentParser()
parser.add_argument("--pop-size",type=int, required=False, dest="ps", default=250, help="Number of diploid individuals")
parser.add_argument("--generations",type=int, required=False, dest="gen", default=50, help="Number of generations")
parser.add_argument("--frequency",type=float, required=False, dest="freq", default=0.5, help="Initial allele frequency")
parser.add_argument("--mcm", action='store_true',default=False,help="Estimate drift distribution with discrete MCM.")
parser.add_argument("--diff", action='store_true',default=False,help="Estimate drift distribution with continuous diffusion model.")
args = parser.parse_args()



ro.r('''

library(ggplot2)


##function for estimating probabilty matrix
drift.prob.nexttgen<- function(proportion=c(), popsize){
  
 
  copies_A<-popsize
  copies_a<-popsize
  total_gametes=copies_A+copies_a
  
  #initialize matrix
   probabilities_matr<- matrix(NA,nrow=total_gametes+1, ncol=total_gametes+1)
   probabilities_matr[1,1]<- 1
   probabilities_matr[nrow(probabilities_matr),ncol(probabilities_matr)]<-1

   
   probabilities=c()
   
  #run binomial distribution to estimate probablities for every possible change
    for (j in 1:(total_gametes-1)){
      p=as.numeric(j/total_gametes)
      q=1-p

        if(p!=0){
      for (y in 0:total_gametes){
        
#        prob<-choose(total_gametes,y) * p^y * q^(total_gametes-y)
	prob<-dbinom(y, total_gametes,p)
        probabilities[y+1]<- abs(prob)
        
      }
        }
	
      if(length(probabilities)>0){
	
      probabilities_matr[j+1,]<-probabilities

      }
   
 
    }
 

  napos2<- which(is.na(probabilities_matr))
  probabilities_matr[napos2]<-0
  return(probabilities_matr)
  
  
    }
  

##run function for estimating probabilities for xxx population size and xxx generations
## make sure that the pop_size is less than or equal to 500, otherwise NAs will be returned

plot_dist<- function(pop_size, generations,freq){
total_allele_counts=2*pop_size
allele_counts=total_allele_counts*freq ##e.g if pop_size=300 (total alleles=600) and denominator=10, then allele_counts=60

allele_freq=allele_counts/total_allele_counts  ##then allele_freq=60/600=0.1
print(paste("Initial frequency of A allele is", allele_freq))

##distribution is polarizes for the rising allele

m<-t(m<-drift.prob.nexttgen(c(rep(0, times=allele_counts), 1,rep(0, times=(total_allele_counts)-allele_counts)), pop_size))
v<-c(rep(0, times=allele_counts), 1,rep(0, times=(total_allele_counts)-allele_counts))


for(i in 1:generations)
{
  v<-c(abs(m%*%v))
}

###plot drift distribution

drift_distr<- matrix(NA, ncol=2, nrow=length(v))
drift_distr[,1]<- v
drift_distr[,2]<- c(0:(length(v)-1))
colnames(drift_distr)<- c("Proportion", "Counts")
drift_distr<- as.data.frame(drift_distr)
return(drift_distr)
}
'''
)







if(args.diff and args.mcm):
	raise ValueError("Don't cheat! Please choose either --mcm or --diff")

elif(args.mcm):
	r_f2=ro.r['plot_dist']
	print("Estimating allele counts")
	res2=r_f2(args.ps,args.gen,args.freq)
#	print(res2)
	print("Start plotting...")
	grdevices = importr('grDevices')

	ro.r('''change_name=function(pop_size, generations,freq){
		name=sprintf("../results/mcm_%sNe_%sfreq_%sgen.png", pop_size,freq,generations)
		return(name)}
	     ''')
	name=ro.r['change_name']
	name=name(args.ps,args.gen,args.freq)
	print("Output figure in:", name)

	grdevices.png(file=name, width=700, height=700)

	gp = ggplot2.ggplot(res2)
	pp=gp+ ggplot2.aes_string(x='Counts', y='Proportion') + ggplot2.geom_bar(stat="identity", color="darkgoldenrod3") + ggplot2.theme_bw()
	pp.plot()
	grdevices.dev_off()
	print("Plot done!")

elif(args.diff): ###references:doi: 10.1093/molbev/msx254 &&  https://doi.org/10.1111/j.1365-294X.2010.04997.x
	p=args.freq
	N=args.ps
	t=args.gen
	newx=[]
	x=np.arange(0,1.001,0.001001001)
	res=[0]*len(x)
	print("Estimating allele counts")
	for i in range(1,101):
		a=p*(1-p)*i*(i+1)*((2*i)+1)
		b=ss.hyp2f1(1-i,i+2,2,p)
		c=math.exp((-i*(i+1)*float(t))/(4*float(N)))
		d=list()
		for xpos in x:
			d.append(ss.hyp2f1(1-i,i+2,2,float(xpos))*a*b*c)
		res=[old + new for old, new in zip(res,d)]

	xcounts=[i*2*N for i in x]
#	print(len(res))
#	print(len(x))
	print("Start plotting...")
	plt.bar(x=x,height=res, width=0.001)
	plt.ylabel('Proportion of population (arbitrary scale)')
	plt.xlabel("Allele Frequency")
	print("Plot done!")
	plt.savefig("../results/diff_{0}Ne_{1}freq_{2}gen.png".format(args.ps,args.freq,args.gen))



elif(not args.diff and not args.mcm):
	 raise ValueError("Please provide either --mcm or --diff, otherwise we will stay here forever!")


else:
	raise ValueError("Something is wrong!!")
