clear all


warning('off','all')

datapath='./../../EmbryoData/Femur_Tibia/';


allpath={
'data/Nuclei_and_CellsE185_S153_m7_distalfemur/',
'data/Nuclei_and_CellsE185_S153_m7_proximaltibia/',
'data/Nuclei_and_CellsE185_S154_m3_distalfemur/',
'data/Nuclei_and_CellsE185_S154_m3_proximaltibia/',
'data/Nuclei_and_CellsE185_S154_m4_distalfemur/',
'data/Nuclei_and_CellsE185_S154_m4_proximaltibia/',
}; 



for gi=1%2:length(allpath)
		path=allpath{gi};
		disp(path)
        s=strsplit(path,'Nuclei_and_CellsE185_');
        outputpath=strcat(datapath,'MakeListNucleiLabelled/',s{2});
        %outputpath_mask=strcat('MakeListClustersMask/',s{2});
        outputpath_unlabelled=strcat(datapath,'MakeListNuclei_unlabelled/',s{2});
        rand_dir_output=strcat(datapath,'RandomSampling/',s{2});
end        


dir2=strcat('RealvsRandomComparison_','/',s{2});
if ~exist([dir2],'dir')
    mkdir([dir2]);
end


realdata=load([outputpath,'centroid_and_surface_nuclei.mat']);
randomdata=load([outputpath_unlabelled,'centroid_and_surface_nuclei.mat']);
%sarah_mask=load([outputpath_mask,'centroid_and_surface_nuclei.mat']);
%sarah_mask_id=load([outputpath,'ClusterMask_NoDuplication.dat']);

[~,LCCreal]=readClusterFile([outputpath,'Cluster_NoDuplication.dat']); 


for realization=1:1
    if mod(realization,10)==0
        realization
    end
    %LCC=readClusterFile(strcat('Realization_DT_E185_nuclei/Random_cluster_',num2str(realization),'.dat'));
    [LCCrandom,allid]=readRandomClusterFile(strcat(rand_dir_output,'Realization/Random_cluster_',num2str(realization),'.dat'));
    %[realization,length(LCC),length(allid)]
    
    randomFeatures=load([rand_dir_output,'AllFeaturesSave',num2str(realization),'.mat']);


%     GVec=averagePC(nuc(allid)); 
%     [Xbin,Ybin_iter, AvgMeanFeatures_iter,AvgStdFeatures_iter]=FindAllFeatures(LCC,rand_dir_output,nuc,c,celvolume,GVec,RZ,HZ,realization,radius);
%     realization_Ybin(:,realization)=Ybin_iter;    
%     for i=1:length(AvgMeanFeatures_iter)
%         realization_AvgMeanFeatures{i}(:,realization)=AvgMeanFeatures_iter{i};
%     end
    
end

realFeatures=load([outputpath,'AllFeaturesSave_.mat']);
sphReal=realFeatures.sphericity; 
rgReal=realFeatures.just_rg; 
p13Real=realFeatures.clusterPC3_by_PC1; 
p23Real=realFeatures.clusterPC3_by_PC2; 



sphRand=randomFeatures.sphericity; 
rgRand=randomFeatures.just_rg; 
p13Rand=randomFeatures.clusterPC3_by_PC1; 
p23Rand=randomFeatures.clusterPC3_by_PC2; 



mycolor={'c.','r.','y.','g.'};

for i=1:length(LCCreal)
   
    idreal=LCCreal{i};
    idrandom=LCCrandom{i};
    
    if length(idreal)>2 
    
    h=figure; 
    set(gcf, 'PaperSize', [10 4]); %7
    set(gcf, 'PaperPosition', [0 0 10 4]);
    subplot(1,2,1)
    for k=1:length(idreal)
        c1=realdata.nuc{idreal(k)};
        ct=mycolor{realdata.unique_tileid(idreal(k),1)};
        plot3(c1(:,1),c1(:,2),c1(:,3),ct,'markersize',1);
        hold on 
    end 
    title(['Real cell=',num2str(length(idreal)) , ',s=', sprintf('%0.2f',sphReal(i)) , ',Rg=', sprintf('%0.1f',rgReal(i)) ,...
        ',PC3/PC1=', sprintf('%0.3f', p13Real(i)) ,',PC3/PC2=', sprintf('%0.3f', p23Real(i))   ],'fontsize',7);
    
    
    subplot(1,2,2)
    for k=1:length(idrandom)
        c1=randomdata.nuc{idrandom(k)};
        
        plot3(c1(:,1),c1(:,2),c1(:,3),'b.','markersize',1);
        hold on 
    end 
    
    title(['Rand cell=',num2str(length(idrandom)) , ',s=', sprintf('%0.2f',sphRand(i)) , ',Rg=', sprintf('%0.1f',rgRand(i)),...
        ',PC3/PC1=', sprintf('%0.3f', p13Rand(i)),',PC3/PC2=', sprintf('%0.3f', p23Rand(i))  ],'fontsize',7);
    saveas(h,[dir2,'cluster','_',num2str(i),'.png']);
    saveas(h,[dir2,'cluster','_',num2str(i),'.fig']);
    close all 
    end
    
end







 function [LCC,allCellId]=readRandomClusterFile(name) 
        fid = fopen(name,'rt');
        tline = fgetl(fid);
        count=1;
        allCellId=[];
         while ischar(tline)
            line= split(tline);
             if length(line)>2
                for j=1:length(line)-1
                    LCC{count}(j)=str2num(line{j});
                    allCellId=[allCellId;LCC{count}(j)];
                end
             end    
%             
%              if length(line)>2
%                 for j=1:length(line)-1
%                     LCC1{count}(j)=str2num(line{j});
%                 end
%              end    
             
            
            tline = fgetl(fid);
            count=count+1;
        end
        fclose(fid);
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

