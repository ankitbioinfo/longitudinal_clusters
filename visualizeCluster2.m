
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
if ~exist([outputVisualize],'dir')
     mkdir([outputVisualize]);
end


for gi=6%2:length(allpath)
		path=allpath{gi};
		disp(path)
        s=strsplit(path,'Nuclei_and_CellsP40_');
        outputpath=strcat('MakeListNucleiLabelled/',s{2});
        
        colmask=load(['MakeListClustersMask/',s{2},'centroid_and_surface_nuclei.mat']);
        %tileid=unique(colmask.unique_tileid(:,3));
        nuclei=load([outputpath,'centroid_and_surface_nuclei.mat']);
        %cloneid=unique(colmask.unique_tileid(:,1));
        
     
               
        savepath=strcat(outputVisualize,s{2});    
        if ~exist([savepath],'dir')
             mkdir([savepath]);
        end

        LCC_mask=load([outputpath,'ClusterMask.dat']);
        
        %many nuclei are duplicated so need to be removed 
        [~,LCC_nuc_duplicated]=readClusterFile([outputpath,'Cluster.dat']);
        clear LCC_nuc
        for i=1:length(LCC_nuc_duplicated)
            id_notgood=LCC_nuc_duplicated{i};
            Repeat_centroids=nuclei.centroid(id_notgood,:);
            [~,firstgood]=unique(Repeat_centroids,'rows');
            id=id_notgood(firstgood);
            Repeat_centroids=nuclei.centroid(id,:);
            Repeat_volume=nuclei.celvolume(id,1);
            clear Repeat_surfaces;
            for j=1:length(id)
                Repeat_surfaces{j}=nuclei.nuc{id(j)};
            end
            
            ia=RemoveBadCells(Repeat_centroids,Repeat_volume,Repeat_surfaces);
            LCC_nuc{i}=id(ia); 
            %[i,length(id),length(ia)]
        end
        
        [LCC_nuc,LCC_mask]=sorted_cluster_from_larger_to_smaller(LCC_nuc,LCC_mask);
       
        fid=fopen([ outputpath,'Cluster_NoDuplication.dat'],'w');
        for i=1:length(LCC_nuc)
            for j=1:length(LCC_nuc{i})
                fprintf(fid,'%d ',LCC_nuc{i}(j));
            end
            fprintf(fid,'\n');
        end
        fclose(fid);
        dlmwrite([outputpath,'ClusterMask_NoDuplication.dat'],LCC_mask);

            clear clen 
            for i=1:length(LCC_nuc)
                id=LCC_nuc{i};
                cent=nuclei.centroid(id,:); 
                t=unique(nuclei.unique_tileid(id,1));
                if length(t)>1
                    for k=1:length(t)
                        temp=find(t(k)==nuclei.unique_tileid(id,1));
                        ankur{i}{k}=id(temp);
                    end
                end
                clen(i,:)=max(cent)-min(cent);
            end

%              [sa,sb]=sort(clen(:,3));
%              dlmwrite('cluster_sorted.dat',[sb,clen(sb,:)],'\t')

            for myid=1:length(LCC_mask)
                id=LCC_nuc{myid};
                if length(id)>5
                    posid=unique(nuclei.unique_tileid(id,3));

                    h=figure;
                    mycolor={'c.','r.','y.'};
                    %yellow replaced by black 
                    ClusterShape=[];
                    mergeVol=0;
                    for i=1:length(id)
                        c=nuclei.nuc{id(i)};
                        %a.unique_tileid(id(i),:)
                        ct=mycolor{nuclei.unique_tileid(id(i),1)};
                        plot3(c(:,1),c(:,2),c(:,3),ct,'markersize',5);
                        hold on 
                        cent=nuclei.centroid(id(i),:);
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
                    title(['# of cells ',num2str(length(id)) , ', section id ',num2str(posid) , ' , cluster height ', num2str(clen(myid,3)) ])

                    %view(6,2); %67
                    %view(-10,-3); %357
                    %view(-96,0) %136
                    %[caz,cel]=view

                    %view(29.4,1.8) %clone 3 cluster 1 
                    %view(14.966,3.02)%clone 3 cluster 2 

                    c=colmask.nuc{LCC_mask(myid)};
                    plot3(c(:,1),c(:,2),c(:,3),'k.','markersize',1);



                    view(-94.64,11.29)


                    saveas(h,[savepath,'cluster',num2str(myid),'.png']);
                    saveas(h,[savepath,'cluster',num2str(myid),'.fig']);
                    close all 
                end
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




function goodcellindex=RemoveBadCells(centroid, volume,surfaces)


        if size(centroid,1)>4
            [~,firstgood]=unique(centroid,'rows');
                %tic
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

            [~,ia]=unique(edges,'rows');
            edges=edges(ia,:);
        else
            edges=[];
            n=size(centroid,1);
            firstgood=1:n;
            for i=1:n 
                for j=i+1:n
                    edges=[edges;[i,j]];
                end
            end
        end    
            
%             edges
%              length(centroid)

             badcell=[];
             for j=1:size(edges,1)
                 dist=pdist(centroid(edges(j,:),:));
                 if dist<10
                         cell1=surfaces{edges(j,1)}; actualVol(1)=volume( edges(j,1)   );      
                         cell2=surfaces{edges(j,2)}; actualVol(2)=volume( edges(j,2)   );  
                         combined=[cell1;cell2];
                         [~,combV]=convhull(combined);
                         [~,comb1]=convhull(cell1);
                         [~,comb2]=convhull(cell2);

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
              %toc
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
