

clear all
allpath={
'data/Nuclei_and_CellsP40_S151_m2_distalfemur/',
'data/Nuclei_and_CellsP40_S151_m2_proximaltibia/',
'data/Nuclei_and_CellsP40_S152_m3_distalfemur/',
'data/Nuclei_and_CellsP40_S152_m3_proximaltibia/',
'data/Nuclei_and_CellsP40_S152_m4_distalfemur/',
'data/Nuclei_and_CellsP40_S152_m4_proximaltibia/',
}; 

nosec={[4,5],[1,2,3],[4,5,6],[1,2,3],[3,4,5],[1,2,3]};
nosec={[4,5],[1,2,3],[4,6],[1,2,3],[3,4,5],[1,2,3]}; %unlabeled 5 is missing 

%allpath={'data/Nuclei_and_Cells_PT_P40_labelled/'}; 
pseudocolor={'unlabeled'}; %'GFP',   
%allpath={'data/Nuclei_and_Cells_PT_P40_unlabelled/'}; pseudocolor={'unlabeled'}; %'GFP',   


 

for gi=3:length(allpath)
		path=allpath{gi};
		disp(path)
        s=strsplit(path,'Nuclei_and_CellsP40_');
        outputpath=strcat('MakeListNuclei_unlabelled/',s{2});    
        if ~exist([outputpath],'dir')
              mkdir([outputpath]);
        end

            
      if exist([outputpath,'centroid_and_surface_nuclei.mat'], 'file') == 0  

            %alignment=load([path,'Alignment_matrix.dat']);
            %vec=alignment(:,1:3);
            %vec=eye(3);
            vec=[0,0,1;0,1,0;1,0,0];


% 
%              [numbers,txt,raw] = xlsread([path,'Tile_coordinates.xlsx']);
%                coordinates = zeros(size(txt,1)-3,5);
%                 for i = 4:size(txt,1),
%                     temp =  char(txt(i,1));
%                     res = strsplit(temp,'_POS');
%                     coordinates(i-3,1) = str2num(char(res(2)));
%                     coordinates(i-3,2:5) = numbers(i-3,:);
%                 end
%                 tile=coordinates(:,2:end);
                
                tile=zeros(20,3);
                tile(:,3)=[0:200:19*200];
                


    clear alltileid
    clear Repeat_surfaces
    %clear Repeat_ellipsoid
    clear Repeat_centroids
    clear Repeat_pixels
    %clear Repeat_triangulation
    clear Repeat_volume
    
%tile=zeros(4,3);    
scount=1;
for clone=1:length(pseudocolor)
for position = nosec{gi}   %coordinates(:,1)'  
    %load(strcat(path,'CreER_masked_nuclei',sprintf('%02d',position),' (Characteristics).mat'));
    sname=s{2}(1:strlength(s{2})-1);
    name=strcat(path,pseudocolor{clone},'_nuclei_',sname,'_sec',sprintf('%d',position), ' (Characteristics).mat');
    load(name);
    %indtemp = find(G.inter.volume_ratio>1);
    %fi = find(coordinates(:,1) == position);
    fi=position; %[fi,tile(fi,:)]
    %for i=1:length(indtemp)
    %    j=indtemp(i);
    for j=1:size(N.coords,1)
            value=G.nuc.centroids(j,:)+ repmat([tile(fi,2), -tile(fi,1),tile(fi,3)],1,1);   
            Repeat_centroids(scount,:)=value*vec;

            value=N.surfaces(j).vertices + repmat([tile(fi,2), -tile(fi,1),tile(fi,3)],1,1);
            Repeat_surfaces{scount,1}=value*vec;
            %Repeat_triangulation{scount,1}=C.surfaces(j).faces;
            %value={G.nuc.ellipsoid_center(j,:), G.nuc.ellipsoid_radii(j,:), G.nuc.ellipsoid_evecs{j}, G.nuc.ellipsoid_v{j}};   
            %Repeat_ellipsoid{scount,1}=value;
            Repeat_ellipsoid(scount,:)=G.nuc.ellipsoid_radii(j,:);
            Repeat_volume(scount,:)=G.nuc.volume(j,:);
            Repeat_pixels(scount,:) = calc_centroids(N.masks(j), N.origins(j,:));
            alltileid(scount,:)=[clone,j,fi];
         
            scount=scount+1;
    end    
end
end


%  ia=RemoveBadCells(Repeat_centroids,Repeat_volume,Repeat_surfaces);
%  [size(Repeat_centroids,1), size(ia,1)]
  ia=1:length(Repeat_centroids);

centroid=Repeat_centroids(ia,:);
nuc=Repeat_surfaces(ia,:);
fitellipsoid=Repeat_ellipsoid(ia,:);
unique_pixel=Repeat_pixels(ia,:);
unique_tileid=alltileid(ia,:);
%faces=Repeat_triangulation(ia,:);
celvolume= Repeat_volume(ia,:);

save([outputpath,'centroid_and_surface_nuclei.mat'],'nuc','centroid','celvolume','unique_pixel','fitellipsoid','unique_tileid','-v7.3');
      else
             load([outputpath,'centroid_and_surface_nuclei.mat'],'centroid');
      end

%         % meanCentroid=mean(centroid);
%         % [meanCentroid, GlobalCenter{gi}{gj}];
%         % [min(centroid(:,3)),max(centroid(:,3))]
%         bonetype=gi;
%         if (bonetype==3)|(bonetype==1)
%                 tempCentroid=[centroid(:,1:2),-centroid(:,3)];
%         else
%                 tempCentroid=[centroid(:,1:2),centroid(:,3)];
%         end
% 
%         newcentroid=tempCentroid-GlobalCenter{gi}{gj};
%         PD_bins=linspace(min(newcentroid(:,3)),max(newcentroid(:,3)),11);
%         
%         if exist([outputpath,'threshold_along_PD_axis.mat'], 'file') == 0  
%         
%             disp('calculate minimum negihbors distance');
%             [min_neigh_dist,clist,edges,remaining_edges,cel_normalizationFactor] = calculate_nuclei_density(newcentroid, [1, 1, 1], 2);
%             individualCellThreshold=min_neigh_dist; 
%             for i=1:10
%                 index=find( (newcentroid(:,3)>PD_bins(i)) &  (newcentroid(:,3)<=PD_bins(i+1)));
%                 median_threshold(i,1)=mean(min_neigh_dist(index));
%                 mean_threshold(i,1)=mean(min_neigh_dist(index));
%             end
% 
%             save([outputpath,'threshold_along_PD_axis.mat'],'median_threshold','mean_threshold','edges','remaining_edges','individualCellThreshold','-v7.3');
%             dlmwrite([outputpath,'threshold_along_PD_axis.dat'],[ median_threshold,   mean_threshold    ],'\t');
%         
%         else
%               load([outputpath,'threshold_along_PD_axis.mat']);
%         end
%         
% 	%[median( individualCellThreshold),mean( individualCellThreshold)]
%         
%         mean_threshold
% 
% % cellsInPZandPHZindex= find((newcentroid(:,3)>-100)&(newcentroid(:,3)<1000));
% % cellsInPZandPHZindex=find(newcentroid(:,3)>PD_bins(3));
% % RZtoHZrange=[min(newcentroid(:,3)),max(newcentroid(:,3))]
% % CellsInPZandPHZrange=[size(centroid,1),length(cellsInPZandPHZindex)]
% disp('I am here')
% 
% fid1=fopen([outputpath,'First_Neighbors.dat'],'w');
% fid2=fopen([outputpath,'Remaining_Neighbors.dat'],'w');
% for i=1:length(centroid)
%     normalized_main=newcentroid(i,:);
%     for j=1:10
%           if( (normalized_main(3)>PD_bins(j)) &  (normalized_main(3)<=PD_bins(j+1)))
%                 median_cutoff=median_threshold(j);
%                 mean_cutoff=mean_threshold(j);
%           end
%     end
%        if   individualCellThreshold(i)<median_cutoff
%       fprintf(fid1,'%d\t%d\n', edges(i,1), edges(i,2));
%        end
%        
%          for j=1:length(remaining_edges{i})
%                 if remaining_edges{i}(j,2)<mean_cutoff
%                     sa=sort([i,remaining_edges{i}(j,1)]);
%                     fprintf(fid2,'%d\t%d\t%f\n', sa(1), sa(2),remaining_edges{i}(j,2));                
%                 end
%          end
%        
%        
%        
% end
% fclose(fid1);
% fclose(fid2);



end


% 
% function [column]=distance_between_neighbor_cell(Cent,ind,cutoff)
%          C=Cent(ind,:);
%          n=size(C,1);
%          c2=1;
%          column=[];
%          for i=1:n
%              for j=i+1:n
%                  ed=pdist(C([i,j],:));  
%                  zd=diff(C([i,j],3));
%                  if ed<=cutoff
%                      column(c2,:)=[ind(i),ind(j),ed,abs(zd)];
%                      c2=c2+1;
%                  end
%              end
%          end
%         
% end



function centroids = calc_centroids(masks, origins)
centroids = nan(length(masks),3);
for i = 1 : length(masks)
    [x,y,z] = ind2sub(size(masks{i}), find(masks{i}));
    centroids(i,:) = (mean([x,y,z],1) + double(origins(i,:) - 1));
end
end






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
     
     dist=pdist(centroid(edges(j,:),:));
     if dist<10   
             cell1=surfaces{edges(j,1)}; actualVol(1)=volume( edges(j,1)   );      
             cell2=surfaces{edges(j,2)}; actualVol(2)=volume( edges(j,2)   );  
             combined=[cell1;cell2];
             [~,comb1]=convhull(cell1);
             [~,comb2]=convhull(cell2);
             [~,combV]=convhull(combined);

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





function LCC=LargestConnectedComponents(edges)
        cellIds=unique(edges(:));
        for j=1:length(cellIds)
            old2new(cellIds(j),1)=j;
            new2old(j,1)=cellIds(j);    
        end
        [length(old2new),length(new2old),length(cellIds)];

        for i=1:size(edges,1)
            for j=1:2 
                newedgename(i,j)= old2new(edges(i,j));
            end
        end
        G=graph(newedgename(:,1),newedgename(:,2));
        bins=conncomp(G);
        % number of connected components 
        nocomp=unique(bins);
        %disp(['# of connected components  ', num2str(length(nocomp))]);
        for i=1:length(nocomp)
            numberOfObjectsInConnectedComponents(i)=sum(nocomp(i)==bins);
        end
        
        
        [sa,sb]=sort(numberOfObjectsInConnectedComponents,'descend');
        
        index=1;
        for i=1:length(sa)
            if sa(i)>1
                LCCIds=find(bins==nocomp(sb(i)));
                LCC{index}=new2old(LCCIds);
                index=index+1;
            end
        end
end




function [neighbor,neighborList,edges1,edges2,convexVolume] = calculate_nuclei_density(N, spacing, delta_spacing, exclude_boundary_points)
%{
Input:
    N - an n*3 matrix with the [x,y,z] coordinates of each of the n points.
    spacing -  1*3 array with the [x,y,z] physical size of the voxels (if N is already in physical units, just use [1,1,1]).
    delta_spacing - a positive real value that specifies the frequency of the sampling of the 3D grid in the output image (the density space / the output argument V).
    exclude_boundary_points - a boolean value that allows the user to choose whether to include or exclude boundary vertices to avoid potential noise in the
    physical boundaries of the sample, default is 'false' (i.e. include boundaries)
Output:
    V - the 3D matrix of estimated densities.
    X, Y, Z - the coordinates of the interpolated space V.
Example:
    N = rand(100, 3);
    spacing = [1 1 1];
    delta_spacing = 0.01;
    exclude_boundary_points = false;
    V = calculate_nuclei_density(N, spacing, delta_spacing, exclude_boundary_points);
    figure;
    imagesc(V(:,:,50));
    axis image;
    colormap jet;
%}
% if the user didn't define 'exclude_boundary_points' we set it to false:
if ~exist('exclude_boundary_points', 'var')
    exclude_boundary_points = false;
end
% calculating the triangulation and the volume of each triangle:
TRI = delaunay(N(:,1), N(:,2), N(:,3));
[~,convexVolume]=convexHull(delaunayTriangulation(N));

clear neighbor
for i = 1 : size(N,1)
    temp=[];
    for j=1:size(TRI,1)
        for k=1:size(TRI,2)
            if TRI(j,k)==i
                temp=[temp,TRI(j,:)];
            end
        end
    end

   
    neighborList{i}=setdiff(unique(temp),i);
    %neighbor(i,1)=length(neighborList{i});
    ids= neighborList{i};
    
    clear dist 
    for k=1:length(ids)
        dist(k)=pdist(N([i,ids(k)],:)); 
    end
   
    [sa,sb]=sort(dist);
    neighbor(i,1)=min(dist);
    
    edges1(i,:) = [sort([i,ids(sb(1))]) sa(1)  ];  
    edges2{i}=[ids(sb(2:end))',sa(2:end)'];
end



%triplot(tri,x,y); for 2d 
%tetramesh(tri,X); 3d 

[size(N,1), length(unique(TRI(:)))];






end
         
