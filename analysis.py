import pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def read_spatial_data(clusterFilename,celltypeFilename):

    df=pd.read_csv(celltypeFilename,sep='\t',header=None)
    data=df.to_numpy()
    spatialcell_unique_clustername=data[:,1]
    spatialcell_unique_clusterid=data[:,0]
    CTname=spatialcell_unique_clustername
    CTid=spatialcell_unique_clusterid

    df=pd.read_csv(clusterFilename)
    louvainFull=df.to_numpy()

    #for i in range(len(CTname)):
    #    print(CTname[i],CTid[i])

    celltype={}
    cellsinCT={}
    index=[]
    for i in range(len(louvainFull)):
        #print(louvainFull[i],louvainFull[i][0])
        clu_id=louvainFull[i][1]
        cel_id=louvainFull[i][0]
        if clu_id in CTid:
            index.append(i)
            #celltype[cel_id]=clu_id
            if clu_id not in cellsinCT:
                cellsinCT[clu_id]=[cel_id]
            else:
                cellsinCT[clu_id].append(cel_id)

    louvain=louvainFull[index,:]
    annotation_spatial_barcode_id= louvain[:,0]
    annotation_spatial_cluster_id= louvain[:,1]



    #print(spatialcell_unique_clustername,spatialcell_unique_clusterid)
    d={}
    for i in range(len(spatialcell_unique_clustername)):
        d[spatialcell_unique_clusterid[i]]=spatialcell_unique_clustername[i]
    annotation_spatial_celltypename=[]
    for i in range(len(annotation_spatial_cluster_id)):
        annotation_spatial_celltypename.append(d[annotation_spatial_cluster_id[i]])
    annotation_spatial_celltypename=np.array(annotation_spatial_celltypename)

    return annotation_spatial_celltypename,annotation_spatial_barcode_id,annotation_spatial_cluster_id,spatialcell_unique_clustername,spatialcell_unique_clusterid



no_of_pc=2
maindir='./'
input_sp_dir='../../data/'

celltypeFilename=input_sp_dir+'NameOfCTmatching.dat'
clusterFilename=input_sp_dir+'sct_leiden_MNN.dat'
#    sct_ad_sp=sc.read_h5ad(fname)
#print(sct_ad_sp)
annotation_spatial_celltypename,annotation_spatial_barcode_id,annotation_spatial_cluster_id,spatialcell_unique_clustername,spatialcell_unique_clusterid=read_spatial_data(clusterFilename,celltypeFilename)


PCA_of_sc_cluster_accordingto_spatial_clusterid=pickle.load(open(maindir+'PCA_of_sc_cluster'+str(no_of_pc)+'.p', 'rb'))
#n=len(input.spatialcell_unique_clustername)


n=len(spatialcell_unique_clustername)
for i in range(n):
    clid=spatialcell_unique_clusterid[i]

    PCA,gene=PCA_of_sc_cluster_accordingto_spatial_clusterid[clid]
    fig,ax=plt.subplots(1,1,figsize=(4,2.5))
    print(PCA.shape)

    ax.plot(PCA[:,0],PCA[:,1],'bo')
    ff=open('fig/'+spatialcell_unique_clustername[i]+'.txt','w')
    sum=0
    for j in range(len(gene)):
        ax.text(PCA[j,0],PCA[j,1],gene[j],fontsize=5)
        ff.write(gene[j]+'\t'+str(PCA[j,0])+'\t'+str(PCA[j,1])+'\n')
        if (np.isnan(PCA[j,0]))|(np.isnan(PCA[j,1])):
            pass
        else:
            sum=sum+PCA[j,0]*PCA[j,1]
    ff.close()

    print(i,sum)

    '''
    fig.tight_layout()
    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')
    #plt.matshow(np.outer(np.sort(model.row_labels_) + 1, np.sort(model.column_labels_) + 1),cmap=plt.cm.Blues)
    #plt.title("Checkerboard structure of rearranged data")

    fig.savefig('fig/corr'+spatialcell_unique_clustername[i]+'.png',dpi=300)
    plt.close('all')
    '''
    #fig.clf()
