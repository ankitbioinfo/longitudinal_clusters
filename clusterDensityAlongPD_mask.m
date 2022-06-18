clear all 
warning('off','all')


allpath={
'data/Nuclei_and_CellsP40_S151_m2_distalfemur/',
'data/Nuclei_and_CellsP40_S151_m2_proximaltibia/',
'data/Nuclei_and_CellsP40_S152_m3_distalfemur/',
'data/Nuclei_and_CellsP40_S152_m3_proximaltibia/',
'data/Nuclei_and_CellsP40_S152_m4_distalfemur/',
'data/Nuclei_and_CellsP40_S152_m4_proximaltibia/',
}; 


RZ_HZ_height={ [-353,511], [-218,635], [-384,1154],[-232,587],[-445,1004],[-227,473]};


for gi=6%2:length(allpath)
		path=allpath{gi};
		disp(path)
        s=strsplit(path,'Nuclei_and_CellsP40_');
        outputpath=strcat('MakeListNucleiLabelled/',s{2});
        outputpath_mask=strcat('MakeListClustersMask/',s{2});
 
end        

load([outputpath,'centroid_and_surface_nuclei.mat']);
sarah_mask=load([outputpath_mask,'centroid_and_surface_nuclei.mat']);
sarah_mask_id=load([outputpath,'ClusterMask_good.dat']);

[~,LCC]=readClusterFile([outputpath,'Cluster_good.dat']); 



bone_mu=mean(centroid);
c=centroid-bone_mu;
% bonetype=1;
% if (bonetype==3)|(bonetype==1)
%        c=[-c(:,1),c(:,2:3)];
%        RotateFactor=[-1 0 0;0 1 0 ;0 0 1];
% else 
%        c=[c(:,1:2),c(:,3)];
% end
RZ=min(c(:,1));
HZ=max(c(:,1));

RotateFactor=eye(3);

% RZ = -299
% HZ = 446
[RZ, HZ]
RZ=RZ_HZ_height{gi}(1);HZ=RZ_HZ_height{gi}(2);
[RZ, HZ]

plot3(c(:,1),c(:,2),c(:,3),'.');
hold on 
axis image
%clear centroid

GVec=averagePC(nuc)  

[vec,val]=eig(cov(c));
mu=mean(c);d = sqrt(diag(val));
factor=3; quiver3(mu(1),mu(2),mu(3),factor*vec(1,1),factor*vec(2,1),factor*vec(3,1),d(1),'r','LineWidth',10);
factor=3; quiver3(mu(1),mu(2),mu(3),factor*vec(1,2),factor*vec(2,2),factor*vec(3,2),d(2),'k','LineWidth',5);
factor=3; quiver3(mu(1),mu(2),mu(3),factor*vec(1,3),factor*vec(2,3),factor*vec(3,3),d(3),'m','LineWidth',2);
 
factor=-3; quiver3(mu(1),mu(2),mu(3),factor*vec(1,1),factor*vec(2,1),factor*vec(3,1),d(1),'r','LineWidth',10);
factor=-3; quiver3(mu(1),mu(2),mu(3),factor*vec(1,2),factor*vec(2,2),factor*vec(3,2),d(2),'k','LineWidth',5);
factor=-3; quiver3(mu(1),mu(2),mu(3),factor*vec(1,3),factor*vec(2,3),factor*vec(3,3),d(3),'m','LineWidth',2);

vec
BoneGrowthAxis=[0; 0 ; 1]
%[0.8709, -0.1407, 0.4708]';


root_dir = '';
addpath(fullfile(root_dir, 'utils'));


dir2=strcat('dataSave2','/');
if ~exist([dir2],'dir')
    mkdir([dir2]);
end



 close all 
 

       

fid=fopen([outputpath,'degree_of_clusters.txt'],'w');

%cluster properties calculate 
 for j=1:length(LCC)
     ClusterShape=[];
     id=LCC{j};
     mergeVol=0;
     for i=1:length(LCC{j})
            %ClusterShape=[  ClusterShape;nuc{LCC{j}(i)}];
            mergeVol=mergeVol+celvolume(LCC{j}(i));
     end
   
     combinedCellVolumeInCluster(j,1)=mergeVol;
     %[K,V]=convhull(ClusterShape); 
     
     
%     plot3(ClusterShape(:,1),ClusterShape(:,2),ClusterShape(:,3),'r.','markersize',5);
%     hold on 
%     [vec,val]=eig(cov(ClusterShape));
%     ovec=vec;
%     oval=val;
%     mu=mean(ClusterShape);
%     d = sqrt(diag(val));
%     hold on;
%     factor=1; quiver3(mu(1),mu(2),mu(3),factor*vec(1,1),factor*vec(2,1),factor*vec(3,1),d(1),'r','LineWidth',10);
%     factor=1; quiver3(mu(1),mu(2),mu(3),factor*vec(1,2),factor*vec(2,2),factor*vec(3,2),d(2),'k','LineWidth',5);
%     factor=1; quiver3(mu(1),mu(2),mu(3),factor*vec(1,3),factor*vec(2,3),factor*vec(3,3),d(3),'m','LineWidth',2);
% 
%     factor=-1; quiver3(mu(1),mu(2),mu(3),factor*vec(1,1),factor*vec(2,1),factor*vec(3,1),d(1),'r','LineWidth',10);
%     factor=-1; quiver3(mu(1),mu(2),mu(3),factor*vec(1,2),factor*vec(2,2),factor*vec(3,2),d(2),'k','LineWidth',5);
%     factor=-1; quiver3(mu(1),mu(2),mu(3),factor*vec(1,3),factor*vec(2,3),factor*vec(3,3),d(3),'m','LineWidth',2);
% 
% 
      mskc=sarah_mask.nuc{sarah_mask_id(j)};
%     plot3(mskc(:,1),mskc(:,2),mskc(:,3),'k.');



     mskvol=sarah_mask.celvolume(sarah_mask_id(j));
     %height=centroid(id,3);
     height=mskc(:,3);
     [ center, radii, evecs, ~, ~ ] = ellipsoid_fit_new( mskc );
     %this step is when fit is not good then use bounding box instead of
     %ellispid fit 
%      if boundingboxVolume < fittedEllipsoidVolume
%          radii = 1/2*(max(ClusterShape)-min(ClusterShape));
%      end
     
      eradii(j,:)=sort(abs(radii)','descend');
      fittedEllipsoidVolume = 4/3*pi*abs(prod(radii)) ;
      boundingboxVolume=  prod(max(mskc)-min(mskc));
      colheight(j,:)=max(height)-min(height);

     just_volume(j,1)=mskvol;
     normalized_volume(j,1)=mskvol/ min (boundingboxVolume, fittedEllipsoidVolume )  ;
     [K,~]=convhull(mskc);  outerSurfacePoints=mskc(unique(K(:)),:);
     just_SA(j,1)=CalculateSurfaceArea(K,mskc);
     normalized_surface_area(j,1)=just_SA(j,1)/ellipsoid_surface_area(radii);
     
     VolumeFraction(j,1)=  mergeVol./mskvol;
    

     sphericity(j,1)= (pi^(1/3))*((6*mskvol).^(2/3)) ./ just_SA(j,1);
     
     clusterSize(j,1)=length(id);
     [clusterPCvector,~,latent]=pca(mskc);
     
     allometric(j,1)=(mskvol^2/3)/just_SA(j,1);
     
     %[ccenter,clusterRadius(j,1)] = sphereFit(ClusterShape);
     %clusterCenter(j,:)=ccenter-mu;
     %clusterRadius(j,1)=   (3*V/(4*pi))^(1/3);
     cluster_mask_fitted_ellipsoid_Rg(j,1)= sqrt(sum(radii.^2)/5);
     
     
     %normalized_rg(j,1)=radiusOfGyration(c,LCC{j}) / cluster_mask_fitted_ellipsoid_Rg(j,1) ;
     just_rg(j,1)=radiusOfGyration2(outerSurfacePoints);
     normalized_rg(j,1)= just_rg(j,1) / cluster_mask_fitted_ellipsoid_Rg(j,1) ;


         
     
     AngleBetweenClusterPC1AndBone_PD(j,1)=180/pi*oangle(clusterPCvector(:,1),BoneGrowthAxis);
     
     clusterPC1(j,1)=latent(1);
     clusterPC2(j,1)=latent(2);
     clusterPC3(j,1)=latent(3);
     
     mainVec=averagePC(nuc(LCC{j}));
     
     clear LNS 
     clear GNS
     clear theta 
     clear biaxial
     
     for i=1:length(LCC{j})
          id=LCC{j}(i);
          [PC,~,~]= pca(nuc{id});
          %oangle(PC(:,1),mainVec(:,1))*180/pi
          LNS(i,1)= 3*(cos(  oangle(PC(:,1),mainVec(:,1))   )^2) -1;
          LNS(i,2)= 3*(cos(  oangle(PC(:,2),mainVec(:,2))  )^2) -1;
          LNS(i,3)= 3*(cos(  oangle(PC(:,3),mainVec(:,3))  )^2)-1;
          
          GNS(i,1)= 3*(cos(  oangle(PC(:,1),GVec(:,1))  )^2) -1;
          GNS(i,2)= 3*(cos(  oangle(PC(:,2),GVec(:,2))  )^2) -1;
          GNS(i,3)= 3*(cos(  oangle(PC(:,3),GVec(:,3))  )^2)-1;     
          
          
          biaxial(i,:)=UniaxialAndBiaxialOOP( PC, mainVec);
          
     end
     
     LOP(j,1)=0.5*(mean(LNS(:,1))); 
     LOP(j,2)=0.5*(mean(LNS(:,2))); 
     LOP(j,3)=0.5*(mean(LNS(:,3))); 
     
     GOP(j,1)=0.5*(mean(GNS(:,1))); 
     GOP(j,2)=0.5*(mean(GNS(:,2))); 
     GOP(j,3)=0.5*(mean(GNS(:,3))); 
     
     LocalBiaxial(j,:)=mean(biaxial,1);

 
     
        [CluAvgDeg(j,1), ColAvgDeg(j,1),largestDistance(j,1),stepsize,deg,edges,AZ_EL_R]=coordination_number(c,LCC{j});
        %coordNumber(j,1)=0; stepsize=0; edges=[]; 
        diameter(j,1)=largestDistance(j,1);
        % this is basically R_g of random walk 
        averageStepSize(j,1) =   sqrt( 1/6*(length(id)*(stepsize^2)));
        AllEdgesOFClusterSaved{j}=AZ_EL_R;
        
        no_ofLoops(j,1)=find_loops_in_cluster(edges);
        
        for i=1:length(deg)
            fprintf(fid,'%d,',deg(i));
        end
        fprintf(fid,'\n');
        
       
        e=normal_modes(c,LCC{j});
        % [length(LCC{j}),nnz(e)]
        highest_mode(j,1)=max(e);
        smallest_mode(j,1)=min(e);
        
        id=LCC{j};
        theta1=[];
        theta2=[];
        theta3=[];
        theta4=[];
        for i=1:size(edges,1)
            flag=[0,0];
            for k=1:length(id)
                if id(k)==edges(i,1)
                    flag(1)=1;
                end
                if id(k)==edges(i,2)
                    flag(2)=1;
                end
            end
            if sum(flag)==2
                 P1= c(edges(i,1),:);  P2= c(edges(i,2),:);
                 P= nuc{edges(i,1)}(1,:);
                 R = cross(P1-P2, P1-P);
                 S = cross(R, P1-P2);
                 alphap=1; betap=1;
                 plane= mean([P1;P2])' + alphap*R' + betap*S';
                 %normCentroidVector=centroidVector/norm(centroidVector);
                 %vec=find_perp(normCentroidVector); vec=vec/norm(vec); [normCentroidVector', vec']
                 %[sum((P1-P2).*S), sum((P1-P2).*R)]
                 
                 
                
                 norm_plane=norm(plane);
                 
                 [PC1,~,~]= pca(nuc{edges(i,1)});  [PC2,~,~]= pca(nuc{edges(i,2)});
                 
                 theta(1)=asin( abs(sum(plane.*PC1(:,1))/norm_plane))*180/pi;
                 theta(2)=asin( abs(sum(plane.*PC2(:,1))/norm_plane))*180/pi;
                 theta1 =[theta1;mean(theta)];
                 
                 theta(1)=asin( abs(sum(plane.*PC1(:,2))/norm_plane))*180/pi;
                 theta(2)=asin( abs(sum(plane.*PC2(:,2))/norm_plane))*180/pi;
                 theta2 =[theta2;mean(theta)];
           
                 theta(1)=asin( abs(sum(plane.*PC1(:,3))/norm_plane))*180/pi;
                 theta(2)=asin( abs(sum(plane.*PC2(:,3))/norm_plane))*180/pi;
                 theta3 = [theta3;mean(theta)];
                 
                 
                 theta4=[theta4;asin( abs(sum(plane.*BoneGrowthAxis)/norm_plane))*180/pi];
                 

            end
        end
        
            
        
        plane_of_division_to_cell(j,:)=[mean(theta1),mean(theta2),mean(theta3)];
        plane_of_division_to_bone(j,1)=mean(theta4);
   
 end
 
 fclose(fid);
 
 fid=fopen([outputpath,'spherical_coordinate.txt'],'w');
 for i=1:length(AllEdgesOFClusterSaved)
     E=AllEdgesOFClusterSaved{i};
     for j=1:size(E,1)
         fprintf(fid,'%d\t%d\t%0.3f\t%0.3f\t%0.3f\t%d\t%d\n',i,E(j,1),E(j,2),E(j,3),E(j,4),E(j,5),E(j,6));
     end
 end
 fclose(fid);

 
clusterPC2_by_PC1=clusterPC2./clusterPC1;
clusterPC3_by_PC1=clusterPC3./clusterPC1;
clusterPC3_by_PC2=clusterPC3./clusterPC2;
% 
% figure
% [sa,sb]=sort(normalized_rg);
% plot(normalized_rg(sb),'b.-')
% hold on 
% plot(clusterRadius(sb),'r.-')
% dlmwrite('ankur.dat',[rg,clusterRadius],'\t')
 
 clusterid=(1:length(eradii))';
 T=table(clusterid,eradii,colheight, CluAvgDeg, ColAvgDeg, clusterPC1,clusterPC2,clusterPC3,normalized_volume , normalized_surface_area,normalized_rg,sphericity, ...
 clusterSize,largestDistance,AngleBetweenClusterPC1AndBone_PD); %angleMeasurement,
 filename = [outputpath,'nuclei_column_stats_data.xlsx'];
 writetable(T,filename,'Sheet','Features','WriteVariableNames',true);
 
 
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 for j=1:length(LCC)
     ClusterShape=[];
     mergeVol=0;
     for i=1:length(LCC{j})
            ClusterShape=[  ClusterShape;nuc{LCC{j}(i)}];
            mergeVol=mergeVol+celvolume(LCC{j}(i));
     end
     warning off;
     [ccenter,clusterRadius1(j,1)] = sphereFit(ClusterShape);
     clusterCenter1(j,:)=(ccenter-bone_mu)*RotateFactor;
     [ccenter2(j,:),tempR] = sphereFit(c(LCC{j},:));
     if tempR>0
        temp_cR(j,1)=min([tempR,clusterRadius1(j,1)]);
     else
        temp_cR(j,1)=clusterRadius1(j,1);
     end
     
 end
 

binsize=51;zone=linspace(RZ,HZ,binsize);Xbin=linspace(0,1,binsize-1);
clusterPositionAndRadius=zeros(length(LCC),2);
for i=1:length(zone)-1
    for j=1:length(LCC)
         cluster= clusterCenter1(j,:);
             if ((cluster(3) >zone(i)) &  (cluster(3) <=zone(i+1)))
                  clusterPositionAndRadius(j,:) = [Xbin(i),  clusterRadius1(j)];
             end
    end
end
dlmwrite([outputpath,'ParameterForRandomModel_basedon_nuc_.dat'], clusterPositionAndRadius,'\t');
% 
% % for i=1:length(zone)-1
% %     for j=1:length(LCC1)
% %          cluster=  ccenter2(j,:);
% %              if ((cluster(1) >zone(i)) &  (cluster(1) <=zone(i+1)))
% %                   clusterPositionAndRadius(j,:) = [Xbin(i),  temp_cR(j)];
% %              end
% %     end
% % end
% % dlmwrite('ParameterForRandomModel_basedon_centroid.dat', clusterPositionAndRadius,'\t');
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  
% 

% 

save([outputpath,'AllFeaturesSave_.mat'],'normalized_volume','normalized_surface_area','clusterSize','sphericity','clusterPC1', 'clusterPC2', 'clusterPC3',...
'clusterPC2_by_PC1','clusterPC3_by_PC1','clusterPC3_by_PC2','just_SA','just_rg',...
    'normalized_rg','VolumeFraction','LOP','GOP','plane_of_division_to_bone', 'no_ofLoops','plane_of_division_to_cell','highest_mode', 'diameter',...
    'CluAvgDeg','smallest_mode','cluster_mask_fitted_ellipsoid_Rg','allometric', 'averageStepSize','LCC','LocalBiaxial','AngleBetweenClusterPC1AndBone_PD','just_volume','colheight');



Features={normalized_volume,just_volume,normalized_surface_area,just_SA,clusterSize,sphericity,clusterPC1, clusterPC2, clusterPC3,clusterPC2_by_PC1,clusterPC3_by_PC1,clusterPC3_by_PC2,...
    normalized_rg,just_rg,VolumeFraction,LOP(:,1),LOP(:,2),LOP(:,3),GOP(:,1),GOP(:,2),GOP(:,3), plane_of_division_to_bone,plane_of_division_to_cell(:,1),...
    plane_of_division_to_cell(:,2),plane_of_division_to_cell(:,3), highest_mode, smallest_mode,cluster_mask_fitted_ellipsoid_Rg, ...,
    LocalBiaxial(:,1), LocalBiaxial(:,2), LocalBiaxial(:,3), LocalBiaxial(:,4), AngleBetweenClusterPC1AndBone_PD,...
    CluAvgDeg,diameter,allometric,averageStepSize,no_ofLoops,colheight};


% savename={'Volume_normalized','Volume','SurfaceArea_normalized','SurfaceArea','clusterSize','clusterSphericity',  ...
%             'clusterPC1','clusterPC2','clusterPC3','clusterPC2_by_PC1','clusterPC3_by_PC1','clusterPC3_by_PC2','radiusOfGyration_normalized','radiusOfGyration',...
%             'volumefraction','LorderPar1','LorderPar2','LorderPar3','GorderPar1','GorderPar2','GorderPar3','plane_of_division_growthAxis',...
%             'plane_of_division_to_PC1', 'plane_of_division_to_PC2', 'plane_of_division_to_PC3','highest_mode', 'smallest_mode','cluster_mask_fitted_ellipsoid_Rg',...
%           'biaxial_S','biaxial_P','biaxial_D','biaxial_C',  'AngleBetweenClusterPC1AndBone_PD' ,'coordination_number','diameter','allometric', 'randomWalktheory','no_ofloops'};
        


% savename{43}='Cluster_Radius';        
% ylabelname={'Cheeger''s Inequality bounds of graphlets','<cluster volume>', '<cluster surface area>','<cluster size>', '<cluster sphericity>',...
%   '<cluster PC1>','<cluster PC2>','<cluster PC3>','<cluster PC2/PC1>','<cluster PC3/PC1>','<cluster PC3/PC2>',...
%  '<conductance>','<r_g>', '<volume fraction>', '<orientational OP1>','<orientational OP2>','<orientational OP3>',...
%  '<orientational OP1>','<orientational OP2>','<orientational OP3>','<E_1>','<E_2>','<E_3>' };
% ylabelname{43}='<cluster radius>';



binsize=21;
zone=linspace(RZ,HZ,binsize);
Xbin=linspace(0,1,binsize-1);
Ybin=zeros(binsize-1,1); 

radius=clusterPositionAndRadius;

for i=1:length(zone)-1
    cellsIdInXInterval=[];
%     for j=1:length(LCC)
%          cluster=mean(c(LCC{j},:));
%              if ((cluster(1) >zone(i)) &  (cluster(1) <=zone(i+1)))
%                  Ybin(i)=Ybin(i)+1;
%                  cellsIdInXInterval=[cellsIdInXInterval,j];
%              end
%     end
    
      for j=1:length(LCC)
           if ((radius(j,1) >Xbin(i)) &  (radius(j,1) <=Xbin(i+1)))
                 Ybin(i)=Ybin(i)+1;
                 cellsIdInXInterval=[cellsIdInXInterval,j];
           end
       end
    
    
    
      for k=1:length(Features)          
                AvgMeanFeatures{k}(i,1)=nanmean(Features{k}(  cellsIdInXInterval,1));
                AvgStdFeatures{k}(i,1)=nanstd(Features{k}(  cellsIdInXInterval,1));
                dataAvgMeanFeatures{k,i}=Features{k}(  cellsIdInXInterval,1);
      end
end


save([outputpath,'SaveVariablesForPlot_.mat'],   'AvgMeanFeatures','AvgStdFeatures','Xbin','Ybin','dataAvgMeanFeatures');


function S=ellipsoid_surface_area(radii)
             p=1.6; a=abs(radii(1));
             b=abs(radii(2)); c=abs(radii(3));
             S=4*pi*(( (a*b)^p + (a*c)^p  + (b*c)^p )/3)^(1/p);
end


function [avgdeg, avgdegExpectedColumn,largestDistance,averageStepSize,deg,E,sph]=coordination_number(centroid,nodeid)
%function [avgdeg,diameter,averageStepSize]=coordination_number(centroids,edges)

%              id=unique(edges(:));
%              for j=1:length(id)
%                  new2old(j)=id(j);
%                  old2new(id(j))=j;
%              end
%                bondDist=[];
%                for i=1:size(edges,1)
%                  newedges(i,1)=old2new(edges(i,1));
%                  newedges(i,2)=old2new(edges(i,2));
%                  bondDist=[bondDist;pdist(centroids(edges(i,:),:))];
%                end
            cent=centroid(nodeid,:);
            n=size(cent,1);
            dist=[];
            myedge=[];
            for i=1:n
                for j=i+1:n
                    dist=[dist,pdist(cent([i,j],:))];
                    myedge=[myedge;[i,j]];
                end
            end
            
            [sa,sb]=sort(dist);
            %[size(myedge), size(dist)]
            
            start=1;
            for i=2:length(sa)
                nodes=unique(myedge(sb(1:i),:));
                % [i,n,length(nodes)]
                if n==length(nodes)
                    start=i;
                    break 
                end
            end
            
          
            
            flag=true;
            while flag
                E=myedge(sb(1:start),:);
                G=graph(E(:,1),E(:,2));
                bins=conncomp(G);
                % number of connected components 
                nocomp=unique(bins);
                clear numberOfObjectsInConnectedComponents
                for i=1:length(nocomp)
                    numberOfObjectsInConnectedComponents(i)=sum(nocomp(i)==bins);
                end

                %[start,numberOfObjectsInConnectedComponents]
                if length(nocomp)==1
                    flag=false;
                    TotalE=start; 
                end
                start=start+1;
            end
               
            
            
            E=myedge(sb(1:TotalE),:);
            averageStepSize=mean(sa(1:TotalE));
            
            %[length(E),length(myedge)]

            largestDistance=sa(TotalE);
            G=graph(E(:,1),E(:,2));
            deg=G.degree();
            avgdegofcompletegraph=length(deg)-1;
            avgdeg=mean(deg);
            avgdegExpectedColumn=2*(n-1)/n;
            %avgdeg= mean(deg)/degreeExpectedColumn; 
            %[length(nodes),length(deg)]
%             
%             d=distances(G);
%             diameter=max(d(:));

%             c=cent-mean(cent);
%             [az,el,r]=cart2sph(c(:,1),c(:,2),c(:,3));
%             sph=[az,el,r];
            
            count=1;
            for i=1:size(E,1)
                c=cent(E(i,:),:);
                c=c-mean(c);
                [az,el,r]=cart2sph(c(:,1),c(:,2),c(:,3));
                flag=[1;1];
                sph(count:count+1,:)=[flag,az,el,r,nodeid(E(i,[1,2]))',nodeid(E(i,[2,1]))'];
                count=count+2;
            end    
            
            for i=1:size(myedge,1)
                c=cent(myedge(i,:),:);
                c=c-mean(c);
                [az,el,r]=cart2sph(c(:,1),c(:,2),c(:,3));
                flag=[0;0];
                sph(count:count+1,:)=[flag,az,el,r,nodeid(myedge(i,[1,2]))',nodeid(myedge(i,[2,1]))'];
                count=count+2;
            end    
                
            
            

end




function e=normal_modes(centroids,nodes)

%          id=unique(edges(:));
%          for j=1:length(id)
%              new2old(j)=id(j);
%              old2new(id(j))=j;
%          end
%          N=length(id);
%          M=zeros(3*N,3*N);
%          for i=1:size(edges,1)
%              r1=centroids(edges(i,1),:);
%              r2=centroids(edges(i,2),:);
%              
%              p1=old2new(edges(i,1));
%              p2=old2new(edges(i,2));
%              
%              q1=3*p1-2; q2=3*p2-2;
%              
%              dist=pdist([r1;r2])^2;
%              
%              M(q1:q1+2,q2:q2+2)=-(1/dist)*   ( (r1-r2)'*(r1-r2) );
%              M(q2:q2+2,q1:q1+2)=-(1/dist)*   ( (r2-r1)'*(r2-r1) );
%          end
%     
         


         id=unique(nodes);
         for j=1:length(id)
             new2old(j)=id(j);
             old2new(id(j))=j;
         end
         N=length(id);
         M=zeros(3*N,3*N);
         for i=1:N
             r1=centroids(id(i),:);
             for j=i+1:N
                     r2=centroids(id(j),:);
             
                     p1=old2new(id(i));
                     p2=old2new(id(j));
             
                     q1=3*p1-2; q2=3*p2-2;

                     dist=pdist([r1;r2])^2;

                     M(q1:q1+2,q2:q2+2)=-(1/dist)*   ( (r1-r2)'*(r1-r2) );
                     M(q2:q2+2,q1:q1+2)=-(1/dist)*   ( (r2-r1)'*(r2-r1) );
             end
         end
         
         
         
         
         %sum(diag(M))
         for i=1:size(M,1)
             M(i,i)=sum(M(i,:));
         end
         
         e=eig(M);
         
end
                



function biaxial=UniaxialAndBiaxialOOP( MoleculePC, LaboratoryPC)
            
         n=MoleculePC(:,1); l=MoleculePC(:,2); m=MoleculePC(:,3); 
         w=LaboratoryPC(:,1);  u=LaboratoryPC(:,2);  v=LaboratoryPC(:,3); 
         
         S = 0.5*  (  (3*(dot(n,w)/(norm(n)*norm(w)) )^2) -1) ; 

         P = 3/2*  (  (dot(n,v)/(norm(n)*norm(v)) )^2 -  (dot(n,u)/(norm(n)*norm(u)) )^2        ) ; 
         
         D = 3/2*  (  (dot(l,w)/(norm(l)*norm(w)) )^2 -  (dot(m,w)/(norm(m)*norm(w)) )^2        ) ;
         
         C=  3/2*  (  (dot(l,v)/(norm(l)*norm(v)) )^2 -  (dot(l,u)/(norm(l)*norm(u)) )^2    +  ...
                      (dot(m,u)/(norm(m)*norm(u)) )^2 -  (dot(m,v)/(norm(m)*norm(v)) )^2                 ) ;
                  
                  
         biaxial=[S,P,D,C];         
end


function rg=radiusOfGyration2(c)
        N= size(c,1);
        center=mean(c,1);
        for i=1:N
            r(i,:)=sum((c(i,:) - center).^2);
        end
        rg=sqrt(sum(r)/N);        
end
      

function rg=radiusOfGyration(centroid,nodes)
        c=centroid(nodes,:);
        N= size(c,1);

        
        center=mean(c,1);
        for i=1:N
            r(i,:)=sum((c(i,:) - center).^2);
        end
        rg=sqrt(sum(r)/N);

%           rg=[]; 
%           for i=1:N
%               for j=1:N
%                   rg=[rg;norm((c(i,:) - c(j,:)))^2  ];
%               end
%           end
%           
%           rg= sqrt(  1/(2*N*N) * sum(rg));
          
          
        
end
            


 function GVec=averagePC(nuc)        
            %length(nuc)
            for i=1:length(nuc)
                [pc,~,~]=pca(nuc{i});
                val1(i,:)=pc(:,1)';
                val2(i,:)=pc(:,2)';
                val3(i,:)=pc(:,3)';
            end

            pc1=pca( [val1;-val1]);
            pc2=pca( [val2;-val2]);
            pc3=pca( [val3;-val3]);

            GVec(:,1)=pc1(:,1);
            GVec(:,2)=pc2(:,1);
            GVec(:,3)=pc3(:,1);
 end




function Area=CalculateSurfaceArea(K,v)
          Area=0;
          for surfind=1:size(K,1)
              pointA=v(K(surfind,1),:);
              pointB=v(K(surfind,2),:);
              pointC=v(K(surfind,3),:);
              TriangleVertex=[pointA; pointB; pointC];
              Ts=pdist(TriangleVertex,'euclidean');
              p=sum(Ts)/2;
              Area=Area+sqrt(p*(p-Ts(1))*(p-Ts(2))*(p-Ts(3)));
          end
          
end


function graphEnergy=calculateEnergy(node,cent)
    
    A=zeros(length(node));
    D=zeros(length(node));
    centroid=cent(node,:);
    for i=1:size(centroid,1)
        for j=i+1:size(centroid,1)
            dist=pdist(centroid([i,j],:));
            A(i,j)= 1/dist;
            A(j,i)= 1/dist;
        end
    end


    for i=1:size(A,1)
        D(i,i)=sum(A(i,:));
    end

    normalizedAdjacency=D^(-1/2)*A*(D^(-1/2));
    L=D-A;
    normalizedLaplacian=D^(-1/2)*L*(D^(-1/2));


    E1=eig(normalizedAdjacency);
    E2=eig(normalizedLaplacian);

    %E2([1:3,end])

    energyAdjacency=sum(abs(E1));

    %totalWeight=sum(sum(A));
    totalWeight=sum(E2);

    energyLaplacian=sum( abs( E2 - (totalWeight/size(A,1))));

    secondSmallestLaplacianEigenvalue=E2(2);

    graphEnergy=[ secondSmallestLaplacianEigenvalue,energyAdjacency,energyLaplacian,min(E1),max(E1)];


end


 function  angle= oangle(u,v)
      angle=atan2(norm(cross(u,v)),dot(u,v));
      if angle>(pi/2)
          angle=pi-angle;
      end
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

 
function no_ofLoops=find_loops_in_cluster(edges)
	G=graph(edges(:,1),edges(:,2));
	t = minspantree(G, 'Type', 'forest');
	%  highlight(h,t)
	nonTreeEdges = setdiff(G.Edges.EndNodes, t.Edges.EndNodes, 'rows');
	cycles_CT = cell(size(nonTreeEdges, 1), 1);
	ankit=[];
	for i = 1 : length(cycles_CT)
	    src = nonTreeEdges(i, 1);
	    tgt = nonTreeEdges(i, 2);
	    cycles_CT{i} = [tgt shortestpath(t, src, tgt)];
	    ankit(i)=length(cycles_CT{i});
	end
	% Christine Tobler's results
	no_ofLoops=length(cycles_CT);
end


                         
function myedges=search_edges(edges,vertices) 
    [~,ia]=unique(edges,'rows');
    edges=edges(ia,:);
	count=1;
    myedges=[];
	for j=1: size(edges,1)
		flag1=0;
		flag2=0;
		for i=1:length(vertices)
			if (edges(j,1)==vertices(i))
				flag1=1;
            end
			if (edges(j,2)==vertices(i))
				flag2=1;
            end
            if (flag1+flag2)==2
                myedges(count,:)=edges(j,:);
                count=count+1;
            end
        end
    end
    
    [~,ia]=unique(myedges,'rows');
    myedges=myedges(ia,:);
    
 end
