
clear all
allpath={
'data/Nuclei_and_CellsP40_S151_m2_distalfemur/',
'data/Nuclei_and_CellsP40_S151_m2_proximaltibia/',
'data/Nuclei_and_CellsP40_S152_m3_distalfemur/',
'data/Nuclei_and_CellsP40_S152_m3_proximaltibia/',
'data/Nuclei_and_CellsP40_S152_m4_distalfemur/',
'data/Nuclei_and_CellsP40_S152_m4_proximaltibia/',
}; 




for gi=1:length(allpath)
		path=allpath{gi};
		disp(path)
        s=strsplit(path,'Nuclei_and_CellsP40_');

        outputpath=strcat('MakeListNucleiLabelled/',s{2});
        
        colmask=load(['MakeListClustersMask/',s{2},'centroid_and_surface_nuclei.mat']);
        tileid=unique(colmask.unique_tileid(:,3));
        nuclei=load([outputpath,'centroid_and_surface_nuclei.mat']);
        cloneid=unique(colmask.unique_tileid(:,1));

        count=1;
        LCC_nuclei=cell(0);
        LCC_mask=[];
        
        for gj=1:length(tileid)
            for gk=1:length(cloneid)
                    [gi,gj,gk]
                    overlap_nuc_id=find(  (nuclei.unique_tileid(:,3)==tileid(gj))&(nuclei.unique_tileid(:,1)==cloneid(gk)) );
                    overlap_mask_id=find((colmask.unique_tileid(:,3)==tileid(gj))&(colmask.unique_tileid(:,1)==cloneid(gk)));
                    
                    [length(overlap_nuc_id),length(overlap_mask_id)]

                    for ii=1:length(overlap_mask_id)
                       
                        i=overlap_mask_id(ii);
                        v1=colmask.nuc{i};
                        [~,comb1]=convhull(v1);
            
%                         plot3(v1(:,1),v1(:,2),v1(:,3),'k.');
%                         hold on 

                 

                                cellincluster=[];
                                for k=1:length(overlap_nuc_id)
                                    j=overlap_nuc_id(k);
                                    v2=nuclei.nuc{j};
                                    %[~,comb2]=convhull(v2);
                                    combined=[v1;v2];
                                    [~,combV]=convhull(combined);
                                    %ankit(k,:)=[comb1,comb2,combV,comb1/combV,j];

                                    %volumeFraction= combV/(comb1+comb2);
                                    vf=comb1/combV;

                                    if vf>0.99
                                        %[gi,gj,gk,count,ii,k]
                                        cellincluster=[cellincluster,j];
                                        %plot3(v2(:,1),v2(:,2),v2(:,3),'b.');
                                    end
                                end
                                if length(cellincluster)>1
                                    LCC_nuclei{count}=cellincluster;
                                    LCC_mask(count,1)=i;
                                    count=count+1;
                                end

                        %end
                    end
            end
        end              

        
    if length(LCC_nuclei)>0    
            [LCC,LCC_mask]=sorted_cluster_from_larger_to_smaller(LCC_nuclei,LCC_mask);
            fid=fopen([ outputpath,'Cluster.dat'],'w');
            for i=1:length(LCC)
                    for j=1:length(LCC{i})
                            fprintf(fid,'%d ',LCC{i}(j));
                    end
                    fprintf(fid,'\n');
            end
            dlmwrite([outputpath,'ClusterMask.dat'],LCC_mask);
            fclose(fid);
    end
       
end



function [LCC,newmask]=sorted_cluster_from_larger_to_smaller(C,mask)
    
    LCC={};
    count=1;
    for i=1:length(C)
         clulen(count,1)=length(C{i});
         count=count+1;
    end
    
    [sa,sb]=sort(clulen,'descend');
    
    for i=1:length(sb)
        LCC{i} = C{ sb(i)};
    end
    
    newmask=mask(sb);

end


