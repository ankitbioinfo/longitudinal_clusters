

clear all
allpath={
'data/Nuclei_and_CellsP40_S151_m2_distalfemur/',
'data/Nuclei_and_CellsP40_S151_m2_proximaltibia/',
'data/Nuclei_and_CellsP40_S152_m3_distalfemur/',
'data/Nuclei_and_CellsP40_S152_m3_proximaltibia/',
'data/Nuclei_and_CellsP40_S152_m4_distalfemur/',
'data/Nuclei_and_CellsP40_S152_m4_proximaltibia/',
}; 


outputVisualize='visualize_nuclei_in_mask_again/';
if ~exist([outputVisualize],'dir')
     mkdir([outputVisualize]);
end




for gi=5   %2:length(allpath)
		path=allpath{gi};
		disp(path)
        s=strsplit(path,'Nuclei_and_CellsP40_');
        outputpath=strcat('MakeListNucleiLabelled/',s{2});
        
        colmask=load(['MakeListClustersMask/',s{2},'centroid_and_surface_nuclei.mat']);
        %tileid=unique(colmask.unique_tileid(:,3));
        a=load([outputpath,'centroid_and_surface_nuclei.mat']);
        %cloneid=unique(colmask.unique_tileid(:,1));
        
end        

savepath=strcat(outputVisualize,s{2});    
        if ~exist([savepath],'dir')
             mkdir([savepath]);
end



 
%[LCC,LCC1]=readClusterFile(['ClusterGreaterThan100Micron1.dat']);
[LCC,LCC1]=readClusterFile([outputpath,'Cluster_good.dat']);
colmaskid=load([outputpath,'ClusterMask_good.dat']);

%col_like_cells=[1496,1526,1561,1584,1614,1657,1696];


for i=1:50%length(LCC)
    id=LCC{i}
    cent=a.centroid(id,:);
    t=unique(a.unique_tileid(id,1));
    if length(t)>1
        for k=1:length(t)
            temp=find(t(k)==a.unique_tileid(id,1));
            ankur{i}{k}=id(temp);
        end
    end
    clen(i,:)=max(cent)-min(cent);
end


for myid= 1:500
id=LCC{myid};
posid=unique(a.unique_tileid(id,3));

h=figure;

mycolor={'c.','r.','y.'};
%yellow replaced by black 
ClusterShape=[];
mergeVol=0;
for i=1:length(id)
    c=a.nuc{id(i)};
    %a.unique_tileid(id(i),:)
    ct=mycolor{a.unique_tileid(id(i),1)};
    flag=0;
    %for j=1:length(col_like_cells)
    %    if id(i)==col_like_cells(j)
    %        flag=1;
    %    end
    %end
    if flag==0
        plot3(c(:,1),c(:,2),c(:,3),ct,'markersize',5);
    else
        plot3(c(:,1),c(:,2),c(:,3),'g.','markersize',5);
    end

    hold on 
    cent=a.centroid(id(i),:);
    %ankit(i,:)=[id(i),cent];
    %ankit(i,:)=max(cent)-min(cent);
    text(cent(1),cent(2),cent(3),num2str(id(i)),'fontsize',8);
    
%     hold on 
%     plot3(cent(:,1),cent(:,2),cent(:,3),'k*')
    
%     Repeat_centroids(i,:)=a.centroid(id(i),:);
%     Repeat_volume(i,:)=a.celvolume(id(i),:);
%     Repeat_surfaces{i}=c;
    
  %   ClusterShape=[  ClusterShape;c];
     %mergeVol=mergeVol+celvolume(LCC{j}(i));
    
end 



%height=max(ankit(:,4))-min(ankit(:,4))
%nocell=length(ankit);
%[max(ankit(:,4)),min(ankit(:,4))]

xlabel('x')
ylabel('y')
zlabel('z')
title(['# of cells ',num2str(length(id)) , ', section id ',num2str(posid)  ])
%title(['# of cells ',num2str(length(id)) , ', section id ',num2str(posid) , ' , cluster height ', num2str(clen(myid,3)) ])
%view(6,2); %67
%view(-10,-3); %357
%view(-96,0) %136
%[caz,cel]=view

%view(29.4,1.8) %clone 3 cluster 1 
%view(14.966,3.02)%clone 3 cluster 2 

c=colmask.nuc{colmaskid(myid)};
plot3(c(:,1),c(:,2),c(:,3),'k.','markersize',1);



view(-94.64,11.29)

saveas(h,[savepath,'cluster',num2str(myid),'.png']);
%saveas(h,[savepath,'cluster',num2str(myid),'.fig']);
close all 

end


 %ia=RemoveBadCells(Repeat_centroids,Repeat_volume,Repeat_surfaces);


function goodcellindex=RemoveBadCells(centroid, volume,surfaces)
    tic
    [~,firstgood]=unique(centroid,'rows');
    %[length(centroid),length(firstgood)]
% spacing=[1, 1, 1];
% delta_spacing=2;
% if the user didn't define 'exclude_boundary_points' we set it to false:
if ~exist('exclude_boundary_points', 'var')
    exclude_boundary_points = false;
end
% calculating the triangulation and the volume of each triangle:
N=centroid(firstgood,:);

TRI = delaunay(N(:,1), N(:,2), N(:,3));
clear neighbor
edges=[];
count=1;
for i = 1 : size(N,1)
    temp=[];
    for j=1:size(TRI,1)
        for k=1:size(TRI,2)
            if TRI(j,k)==i
                temp=[temp,TRI(j,:)];
            end
        end
    end

   
    neighborList=setdiff(unique(temp),i);
    %neighbor(i,1)=length(neighborList{i});
    ids= neighborList;
    
    for j=1:length(ids)
        a=min(ids(j),i);
        b=max(ids(j),i);
        edges(count,:)=[firstgood(a),firstgood(b)];
        count=count+1;
    end
 
end

%size(edges)
[~,ia]=unique(edges,'rows');
edges=edges(ia,:);
%size(edges)

 badcell=[];
 for j=1:size(edges,1)
     edges(j,:)
     dist=pdist(centroid(edges(j,:),:));
     if dist<10
             cell1=surfaces{edges(j,1)}; actualVol(1)=volume( edges(j,1)   );      
             cell2=surfaces{edges(j,2)}; actualVol(2)=volume( edges(j,2)   );  
             combined=[cell1;cell2];
             [~,combV]=convhull(combined);
             [~,comb1]=convhull(cell1);
             [~,comb2]=convhull(cell2)

             volumeFraction= combV/(comb1+comb2);
             if volumeFraction<1
                        [sa,sb]=min(actualVol);
                        if sb==1
                            badid=edges(j,1);
                        end
                        if sb==2
                            badid=edges(j,2);
                        end
                        badcell=[badcell;badid];
             end
     end
 end    

  goodcellindex=setdiff(firstgood,badcell);
  toc
  %length(goodcellindex)  
  
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
