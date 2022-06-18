clear all 


allpath={
'data/Nuclei_and_CellsP40_S151_m2_distalfemur/',
'data/Nuclei_and_CellsP40_S151_m2_proximaltibia/',
'data/Nuclei_and_CellsP40_S152_m3_distalfemur/',
'data/Nuclei_and_CellsP40_S152_m3_proximaltibia/',
'data/Nuclei_and_CellsP40_S152_m4_distalfemur/',
'data/Nuclei_and_CellsP40_S152_m4_proximaltibia/',
}; 

outputVisualize='visualize_nuclei_in_mask/';

for gi=1:6
    path=allpath{gi};
		disp(path)
        s=strsplit(path,'Nuclei_and_CellsP40_');
        outputpath=strcat('MakeListNucleiLabelled/',s{2});
        outputpath_mask=strcat('MakeListClustersMask/',s{2});
        
        load([outputpath,'centroid_and_surface_nuclei.mat']);
        sarah_mask=load([outputpath_mask,'centroid_and_surface_nuclei.mat']);
        sarah_mask_id=load([outputpath,'ClusterMask_NoDuplication.dat']);
        [~,LCC]=readClusterFile([outputpath,'Cluster_NoDuplication.dat']); 
        
        %sarah_mask_id=load([outputpath,'ClusterMask.dat']);
        %[~,LCC]=readClusterFile([outputpath,'Cluster.dat']); 
        
         for j=1:length(LCC)
            mskc=sarah_mask.nuc{sarah_mask_id(j)};
            %plot3(mskc(:,1),mskc(:,2),mskc(:,3),'k.');
            %mskvol(j,1)=sarah_mask.celvolume(sarah_mask_id(j));
            [K,V]=convhull(mskc); 
            mskvol(j,1)=V;
            height=mskc(:,3);
            %colheight(j,1)=max(height)-min(height);
            colheight(j,1)= length(LCC{j});
         end
         
         [sa,sb]=sort(mskvol,'descend');
          %savepath=strcat(outputVisualize,s{2});    
         dlmwrite( strcat(outputVisualize,'maskVol',num2str(gi),'.dat'),[sb,sa, colheight(sb)],'\t');
         clear mskvol
         clear colheight
end



 
 function [LCC,LCC1]=readClusterFile(name) 
        fid = fopen(name,'rt');
        tline = fgetl(fid);
        count=1;
         while ischar(tline)
            line= split(tline,' ');
             if length(line)>3
                for j=1:length(line)-1
                    LCC{count}(j)=str2num(line{j});
                end
             end    
            
             if length(line)>2
                for j=1:length(line)-1
                    LCC1{count}(j)=str2num(line{j});
                end
             end    
             
            
            tline = fgetl(fid);
            count=count+1;
        end
        fclose(fid);
 end
