clear all 


allpath={
'data/Nuclei_and_CellsP40_S151_m2_distalfemur/',
'data/Nuclei_and_CellsP40_S151_m2_proximaltibia/',
'data/Nuclei_and_CellsP40_S152_m3_distalfemur/',
'data/Nuclei_and_CellsP40_S152_m3_proximaltibia/',
'data/Nuclei_and_CellsP40_S152_m4_distalfemur/',
'data/Nuclei_and_CellsP40_S152_m4_proximaltibia/',
}; 


outputVisualize='cartoon_figure_for_EL_AZ/';
if ~exist([outputVisualize],'dir')
     mkdir([outputVisualize]);
end



for gi=6   %2:length(allpath)
		path=allpath{gi};
		disp(path)
        s=strsplit(path,'Nuclei_and_CellsP40_');
        outputpath=strcat('MakeListNucleiLabelled/',s{2});
        
        colmask=load(['MakeListClustersMask/',s{2},'centroid_and_surface_nuclei.mat']);
        %tileid=unique(colmask.unique_tileid(:,3));
        a=load([outputpath,'centroid_and_surface_nuclei.mat']);
        %cloneid=unique(colmask.unique_tileid(:,1));
        
end  

[LCC,LCC1]=readClusterFile([outputpath,'Cluster_good.dat']);
colmaskid=load([outputpath,'ClusterMask_good.dat']);


bone_mu=mean(a.centroid);
gcent=a.centroid-bone_mu;



for myid= 1011:1269
id=LCC1{myid};
posid=unique(a.unique_tileid(id,3));

h=figure;

mycolor={'c.','r.','y.'};
mycolor={'b.','b.','b.'};
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
    
    
end
    
   [AZ_EL_R]=coordination_number(gcent,id);
   
   alphabet='ABCDEFGHIJ';
   
   cu=a.centroid(id,:); %cu=cu-mean(cu); 
   cube=[cu(1,1),cu(1,2),cu(1,3)
        cu(2,1),cu(1,2),cu(1,3)
        cu(2,1),cu(2,2),cu(1,3)
        cu(1,1),cu(2,2),cu(1,3)
        cu(1,1),cu(1,2),cu(1,3)];
   plot3(cube(:,1),cube(:,2),cube(:,3),'k-');
   A=cube(1,:);    B=cube(2,:); C=cube(3,:);
   AB=B-A; AC=C-A; P=cube(4,:);
   
   for i=1:4
       text(cube(i,1),cube(i,2),cube(i,3),alphabet(i),'fontsize',12);
   end
   cube=[cu(1,1),cu(1,2),cu(2,3)
        cu(2,1),cu(1,2),cu(2,3)
        cu(2,1),cu(2,2),cu(2,3)
        cu(1,1),cu(2,2),cu(2,3)
        cu(1,1),cu(1,2),cu(2,3)];
    plot3(cube(:,1),cube(:,2),cube(:,3),'k-');
   for i=1:4
       text(cube(i,1),cube(i,2),cube(i,3),alphabet(i+4),'fontsize',12);
   end
    
    
    cube=[cu(1,1),cu(1,2),cu(1,3); cu(1,1),cu(1,2),cu(2,3)];plot3(cube(:,1),cube(:,2),cube(:,3),'k-');
    cube=[cu(2,1),cu(1,2),cu(1,3); cu(2,1),cu(1,2),cu(2,3)];plot3(cube(:,1),cube(:,2),cube(:,3),'k-');
    cube=[cu(2,1),cu(2,2),cu(1,3); cu(2,1),cu(2,2),cu(2,3)];plot3(cube(:,1),cube(:,2),cube(:,3),'k-');
    cube=[cu(1,1),cu(2,2),cu(1,3); cu(1,1),cu(2,2),cu(2,3)];plot3(cube(:,1),cube(:,2),cube(:,3),'k-');
    cube=[cu(1,1),cu(1,2),cu(1,3); cu(2,1),cu(2,2),cu(2,3)];plot3(cube(:,1),cube(:,2),cube(:,3),'k-','linewidth',2);

    pointP=cu(2,:);
    n = cross(AB, AC) ; 
    n = n / norm(n) ;
    % project onto the plane
    C_proj = pointP - dot(pointP - P, n) * n; 
    
    cube=[A;C_proj]; plot3(cube(:,1),cube(:,2),cube(:,3),'k:','linewidth',2);
    
    x_axis=[1,0,0];  GA=diff(cu); proj=C_proj-A;
    AZ= oangle(x_axis',(proj)');
    EL = oangle(GA', proj');
    [AZ, EL]
    
   
   
   
   %2 is AZ and 3 is EL 
    
    %index=find(AZ_EL_R(:,3)>0);
    %mydata=AZ_EL_R(index,:);
    
    for i=1:size(AZ_EL_R,1)
        cid=AZ_EL_R(i,5);
        cent=a.centroid(cid,:);
        %ankit(i,:)=[id(i),cent];
        %ankit(i,:)=max(cent)-min(cent);
        ss=strcat('  ',num2str(cid), ',[AZ=',sprintf('%0.1f',AZ_EL_R(i,2)) , ', EL=',sprintf('%0.1f',AZ_EL_R(i,3)), ', r=', sprintf('%0.1f',AZ_EL_R(i,4)) , ']'  );
        text(cent(1),cent(2),cent(3),ss,'fontsize',14);
    end
    
    
    
    
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

    %c=colmask.nuc{colmaskid(myid)};
    %plot3(c(:,1),c(:,2),c(:,3),'k.','markersize',1);


    view(-94.64,11.29)

    saveas(h,[outputVisualize,'cluster',num2str(myid),'.png']);
    saveas(h,[outputVisualize,'cluster',num2str(myid),'.fig']);
    close all 

end



function [sph]=coordination_number(centroid,nodeid)
                    factor=180/pi;
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
                sph(count:count+1,:)=[flag,az*factor,el*factor,r,nodeid(E(i,[1,2]))',nodeid(E(i,[2,1]))'];
                count=count+2;
            end    
            
%             for i=1:size(myedge,1)
%                 c=cent(myedge(i,:),:);
%                 c=c-mean(c);
%                 [az,el,r]=cart2sph(c(:,1),c(:,2),c(:,3));
%                 flag=[0;0];
%                 sph(count:count+1,:)=[flag,az,el,r,nodeid(myedge(i,[1,2]))',nodeid(myedge(i,[2,1]))'];
%                 count=count+2;
%             end    
                
            
            

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