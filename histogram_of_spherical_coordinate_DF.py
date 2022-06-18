
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
from scipy.optimize import curve_fit
from scipy.stats import genextreme as gev


def bimodal(x,mu1,mu2,s1,s2,p1):
    disA=np.exp(-0.5*(x-mu1)**2 / s1**2.) / (s1*np.sqrt(2.0 * np.pi))
    disB=np.exp(-0.5*(x-mu2)**2 / s2**2.) / (s2*np.sqrt(2.0 * np.pi))
    return p1*disA+(1-p1)*disB

def gumbel(x,mu,beta):
    return np.exp(-np.exp(-(x-mu)/beta))


def lognormal(x,mu,sigma):
	return (1.0/(x*sigma*np.sqrt(2*np.pi)))  *np.exp(- (np.log(x)-mu)**2 / (2*(sigma**2)) )


def exponential(x, a, b, c):
     return a*np.exp(-b*x) + c


def plot_everything(datapath,ax,mycolor,legendlabel):

	f=open(datapath+'degree_of_clusters.txt')

	deg=[]
	for line in f:
		l=line.split(',')
		for i in range(len(l)-1):
			deg.append(int(i))	



	data=pd.read_csv(datapath+'spherical_coordinate.txt',sep='\t',header=None)
	#data=pd.read_csv('Embryo/spherical_coordinate12.txt',sep='\t',header=None)

	data=data.to_numpy()
	print(data.shape)

	factor=180/np.pi



	AZ=data[:,1]*factor
	#legend= ax.legend(loc='best',ncol=1, borderaxespad=0., prop={"size":6},fancybox=True, shadow=True)
	ax[0,0].hist(AZ,bins=50,density=True,histtype='step', color=mycolor, facecolor=mycolor,alpha=1,label=legendlabel)

	ax[0,0].set_xlim([-180,180.1])
	ax[0,0].set_xlabel('Azimuthal')
	ax[0,0].set_ylabel('P(Azimuthal)')
	start, end = ax[0,0].get_xlim()
	stepsize=60
	ax[0,0].xaxis.set_ticks(np.arange(np.floor(start), end, stepsize))




	x=data[:,3]
	n,bins,_=ax[0,1].hist(x,bins=50,density=True,histtype='step', color=mycolor,facecolor=mycolor,alpha=1)

	'''
	hist=ax[0,1].hist(x,bins=50,density=True,facecolor='g',alpha=0.75)
	x = [hist[0], 0.5*(hist[1][1:]+hist[1][:-1])]
	xdata=x[1]
	ydata=x[0]
	#plt.plot(xdata,ydata,'k-')
	popt, pcov = curve_fit(gumbel, xdata, ydata)
	#print('point',popt)
	ax[0,1].plot(xdata, gumbel(xdata, *popt), 'k-', label='fit: \mu=%5.3f, \beta=%5.3f' % tuple(popt))
	'''
	shape,loc,scale = stats.lognorm.fit(x)
	print(shape,loc,scale)
	mu=np.log(scale)
	sigma=shape
	best_fit_line = stats.lognorm.pdf(bins,shape, loc=loc, scale=scale)
	#ax[0,1].set_title(r'$\mu$='+'%0.2f'%mu+', $\sigma$='+'%0.2f'%sigma+', loc='+'%0.2f'%loc)
	legend2=r':$\mu$='+'%0.2f'%mu+', $\sigma$='+'%0.2f'%sigma+', loc='+'%0.2f'%loc
	ax[0,1].plot(bins, best_fit_line,mycolor+'-',label=legendlabel+legend2)

	#ax[0,1].plot(bins,lognormal(bins,mu,shape),'ro')
	#popt, pcov = curve_fit(lognormal, xdata, ydata)
	#print(popt,pcov)

	ax[0,1].set_xlim([0,25])
	ax[0,1].set_xlabel('r')
	ax[0,1].set_ylabel('P(r)')

	EL=data[:,2]*factor
	ax[1,0].hist(EL,bins=50,density=True,histtype='step',color=mycolor, facecolor=mycolor,alpha=1,label=legendlabel)
	ax[1,0].set_xlim([-90,90.1])
	ax[1,0].set_xlabel('Elevation')
	ax[1,0].set_ylabel('P(Elevation)')
	start, end = ax[1,0].get_xlim()
	stepsize=30
	ax[1,0].xaxis.set_ticks(np.arange(np.floor(start), end, stepsize))


	'''
	#nodes=data[:,[4,5,6]]
	count=0
	edges=[]
	for i in range(len(AZ)):
		if abs(EL[i])>60:
			#if ((abs(AZ[i]) < 30) | (abs(AZ[i]) > 150)):  
				count=count+1
				edges.append([str(int(data[i,6]))+'_'+str(int(data[i,4])), str(int(data[i,6]))+'_'+str(int(data[i,5]))])
				
	print(edges)			

	print('AZ', len(AZ), 'EL', len(EL), 'GT', count, 100*count/len(AZ))
	'''

	bins=np.arange(min(deg), max(deg)+1)-0.5
	hist=ax[1,1].hist(deg,bins=bins,density=True,histtype='step', color=mycolor,facecolor=mycolor,alpha=1)
	x = [hist[0], 0.5*(hist[1][1:]+hist[1][:-1])]
	xdata=x[1]
	ydata=x[0]
	popt, pcov = curve_fit(exponential, xdata, ydata)
	legend4=':%5.3f exp(-%5.3fx) + %.2f' % tuple(popt)
	ax[1,1].plot(xdata, exponential(xdata, *popt), mycolor+'-',label=legendlabel+legend4)

	#ax[1,1].set_title('%5.3f exp(-%5.3fx) + %.2f' % tuple(popt)) 
	ax[1,1].set_xlabel('deg')
	ax[1,1].set_ylabel('P(deg)')
	ax[1,1].set_xlim([0,40])
	
	#return legend2,lengend4


def main():
	datapath1='plothist_embryo/DF/'
	datapath2='plothist_postnatal/DF/'


	fig,ax=plt.subplots(2,2,figsize=(8,5))
	
	plot_everything(datapath1,ax,'r','E18.5 DF')
	plot_everything(datapath2,ax,'b','P40 DF')
	ax[0,0].legend(prop={'size': 7},loc='lower center')
	ax[0,1].legend(prop={'size': 6},loc='upper right')
	ax[1,0].legend(prop={'size': 7},loc='lower center')
	ax[1,1].legend(prop={'size': 6},loc='upper right')

	fig.tight_layout()
	#fig.savefig('Embryo_automated_cluster_cutoff_12.png',bbox_inches='tight',dpi=300)
	fig.savefig('histogram_elevation_cluster_DF.png',bbox_inches='tight',dpi=300)
	fig.clf()
	
main()	
