
clear all 

%datafile='plothist_embryo/PT/';
datafile='plothist_embryo/DF/';

data=load([datafile,'spherical_coordinate.txt']);
data_EL_conn_comp=load([datafile,'spherical_coordinate_not_in_connected_component.txt']);


factor=180/pi;
EL=data(:,3)*factor;

cutoff=60;

count=0;
column_like=[];
for i=1:length(EL)
    if abs(EL(i))>cutoff
        count=count+1;
        column_like(count,:)=data(i,[5,6]);
    end
end

temp=data(:,[5,6]);
total_no_of_cells=length(unique(temp(:)))
total_no_of_cells_in_col=length(unique(column_like(:)))

%doublet_col=[100*count/length(EL)]

EL_doublet_unique=length(unique(column_like,'rows'))
all_doublet_unique=length(  unique(data(:,[5,6]),'rows' )  )
total_no_of_doublets=100*EL_doublet_unique/all_doublet_unique

cells_in_doublet=100*total_no_of_cells_in_col/total_no_of_cells

LCC=LargestConnectedComponents(column_like);

length(LCC)



columnlike=zeros(1,6);
clusterlike=zeros(1,6);

for i=1:length(LCC)
    if length(LCC{i})==2
        columnlike(2)=columnlike(2)+1;
    end
    
    if length(LCC{i})==3
        
        if find_real_columns_or_not(LCC{i},data_EL_conn_comp,cutoff)
                columnlike(3)=columnlike(3)+1;
        else
                clusterlike(3)=clusterlike(3)+1;
        end
    end
    
    if length(LCC{i})==4
        if find_real_columns_or_not(LCC{i},data_EL_conn_comp,cutoff)
            columnlike(4)=columnlike(4)+1;
        else
            clusterlike(4)=clusterlike(4)+1;

        end
    end
    
    if length(LCC{i})>5
        if find_real_columns_or_not(LCC{i},data_EL_conn_comp,cutoff)
            columnlike(6)=columnlike(6)+1;
            
        else
            clusterlike(6)=clusterlike(6)+1;
        end
    end
    
     if length(LCC{i})==5
          if find_real_columns_or_not(LCC{i},data_EL_conn_comp,cutoff)
                columnlike(5)=columnlike(5)+1;
          else
                clusterlike(5)=clusterlike(5)+1;

          end
     end
    
  
    
  
    
end
        
collen=100*columnlike/length(LCC);
cullen=100*clusterlike/length(LCC);

columnlike
clusterlike



function flag=find_real_columns_or_not(nodes,data,cutoff)
        factor=180/pi;
        EL=data(:,3)*factor;
        %nodes
        count=0;
        total=0;
        for i=1:length(nodes)
            for j=i+1:length(nodes)
                index=find(  (data(:,5)==nodes(i))&(data(:,6)==nodes(j)));
                total=total+1;
                if abs(EL(index,:))>cutoff
                    count=count+1;
                end
            end
        end
        
        if total==count
            flag=true;
        else
            flag=false;
        end
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
