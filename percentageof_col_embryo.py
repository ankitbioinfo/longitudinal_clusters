import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
plt.rcParams.update({'font.size': 12})


myxlabel=['2','3','4','5','>5']
tname=['P40 PT (EL >70)','P40 DF (EL >70)','P40 PT (EL >60)','P40 DF (EL >60)']

fig, ax = plt.subplots(2,2, figsize=(10, 7))	
#print('embryo',readfiles(name1,'NA'),readfiles(name11,'NA'))

total=373
a=[339,12,4,0,0]
b=[0,14,4,0,0]
print('1',sum(a)+sum(b),total)
a=np.array(a)/total
b=np.array(b)/total
ratio={'Cluster-like':b,'Column-like':a}
df=pd.DataFrame(ratio,index=myxlabel)
df.plot(ax=ax[0,0],kind='bar',stacked=True)
ax[0][0].set_title(tname[0])
#ax[0][0].text(0,1,c[0])
#ax[0][0].text(1,1,c[1])
ax[0][0].set_ylabel('frequency')
ax[0][0].legend(loc='lower left',fontsize=8, bbox_to_anchor=(0.7,1))




total=317
a=[290,9,1,0,0]
b=[0,14,3,0,0]
print('2',sum(a)+sum(b),total)
a=np.array(a)/total
b=np.array(b)/total
ratio={'Cluster-like':b,'Column-like':a}
df1=pd.DataFrame(ratio,index=myxlabel)
df1.plot(ax=ax[0,1],kind='bar',stacked=True)
ax[0][1].set_title(tname[1])
#ax[0][1].text(0,1,c[0])
#ax[0][1].text(1,1,c[1])
ax[0][1].legend(loc='lower left',fontsize=8, bbox_to_anchor=(0.7,1))



total=627
a=[500,46,6,0,0]
b=[0,34,23,5,13]
print('3',sum(a)+sum(b),total)
a=np.array(a)/total
b=np.array(b)/total
ratio={'Cluster-like':b,'Column-like':a}
df2=pd.DataFrame(ratio,index=myxlabel)
df2.plot(ax=ax[1,0],kind='bar',stacked=True)
ax[1][0].set_title(tname[2])
#ax[1][0].text(0,1,c[0])
#ax[1][0].text(1,1,c[1])
ax[1][0].set_ylabel('frequency')
ax[1][0].set_xlabel('Number of cells in columns/clusters')
ax[1][0].legend(loc='lower left',fontsize=8, bbox_to_anchor=(0.7,1))


total=533
a=[461,24,0,0,0]
b=[0,22,15,5,6]
print('4',sum(a)+sum(b),total)
a=np.array(a)/total
b=np.array(b)/total

ratio={'Cluster-like':b,'Column-like':a}
df3=pd.DataFrame(ratio,index=myxlabel)
df3.plot(ax=ax[1,1],kind='bar',stacked=True)
ax[1][1].set_title(tname[3])
#ax[1][1].text(0,1,c[0])
ax[1][1].set_xlabel('Number of cells in columns/clusters')
ax[1][1].legend(loc='lower left',fontsize=8, bbox_to_anchor=(0.7,1))


fig.tight_layout()
fig.savefig('Ideal_columns_percentage_embryo.png',bbox_inches='tight',dpi=300)

