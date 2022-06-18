clear all 
 
 
% 1-5  CheegerInequality,volume,surface_area,clusterSize,sphericity,
% 6-10 clusterPC1, clusterPC2, clusterPC3,clusterPC2_by_PC1,clusterPC3_by_PC1,
% 11-14 clusterPC3_by_PC2,cycle,diameter,Res3};


% 1: index 
% 2: diameter 
% 3: transitivity 
% 4: clustering coefficient 1
% 5: clustering coefficient 2
% 6: average degree 1
% 7: average degree 2
% 8: efficiency global 
% 9: efficiency local 
% 10: wiener 1 
% 11: wiener 2 
% 12: estrada_index
% 13: subgraph_centrality
% 14: node_connectivity 
% 15: edge_connectivity 
% 16: current_flow_betweenness_centrality_node_1 
% 17: current_flow_betweenness_centrality_node_2
% 18: current_flow_betweenness_centrality_edge_1 
% 19: current_flow_betweenness_centrality_edge_2
% 20: resistance_distance1
% 21: resistance_distance2
% 22: information_centrality1
% 23: information_centrality2
% 24: average degree connectivity 1 
% 25: average degree connectiivty 2 
% 26: robustness 



savename={'Volume_normalized','Volume','SurfaceArea_normalized','SurfaceArea','clusterSize','clusterSphericity',  ...
            'clusterPC1','clusterPC2','clusterPC3','clusterPC2_by_PC1','clusterPC3_by_PC1','clusterPC3_by_PC2','radiusOfGyration_normalized','radiusOfGyration',...
            'volumefraction','LorderPar1','LorderPar2','LorderPar3','GorderPar1','GorderPar2','GorderPar3','plane_of_division_growthAxis',...
            'plane_of_division_to_PC1', 'plane_of_division_to_PC2', 'plane_of_division_to_PC3','highest_mode', 'smallest_mode','cluster_mask_fitted_ellipsoid_Rg',...
          'biaxial_S','biaxial_P','biaxial_D','biaxial_C',  'AngleBetweenClusterPC1AndBone_PD' ,'coordination_number','diameter','allometric', 'randomWalktheory','no_ofloops','columnHeight'};

ylabelname={'$\langle\overline{cluster\; volume}\rangle$', '<cluster volume>',  '$\langle\overline{cluster\; surface\; area}\rangle$',  '<cluster surface area>',  '<cluster size>', '<cluster sphericity>',...
  '<cluster PC1>','<cluster PC2>','<cluster PC3>','<cluster PC2/PC1>','<cluster PC3/PC1>','<cluster PC3/PC2>', '$\langle\overline{R_g}\rangle$','<R_g>',...
   '<volume fraction>', '<Local OOP 1>','<Local OOP 2>','Local <OOP 3>',...
   '<Global OOP 1>','<Global OOP 2>','<Global OOP 3>',  '<angle(plane, P-D axis)>','<angle(plane,PC1)>',...
            '<angle(plane, PC2)>', '<angle(plane,PC3)>', '<largest eigenvalue of Hessian>', '<smallest eigenvalue of Hessian>',...
             '<R_g>', '< OOP S>', '< OOP P>', '< OOP D>', '< OOP C>', '<\alpha>','$\langle{avg\;deg}\rangle$','<diameter>','Vol^{2/3}/SA','<R_g>','<# of cycles>','<column height>'};
        
         
latexSymbol=[1,3,13,34];         
        

dir2=strcat('ClusterFeatures_profiles_PT_E185_P40','/');
if ~exist([dir2],'dir')
    mkdir([dir2]);
end


allpathEmbryo={
'data/Nuclei_and_CellsE185_S153_m7_distalfemur/',
'data/Nuclei_and_CellsE185_S153_m7_proximaltibia/',
'data/Nuclei_and_CellsE185_S154_m3_distalfemur/',
'data/Nuclei_and_CellsE185_S154_m3_proximaltibia/',
'data/Nuclei_and_CellsE185_S154_m4_distalfemur/',
'data/Nuclei_and_CellsE185_S154_m4_proximaltibia/',
}; 

allpathPostnatal={
'data/Nuclei_and_CellsP40_S151_m2_distalfemur/',
'data/Nuclei_and_CellsP40_S151_m2_proximaltibia/',
'data/Nuclei_and_CellsP40_S152_m3_distalfemur/',
'data/Nuclei_and_CellsP40_S152_m3_proximaltibia/',
'data/Nuclei_and_CellsP40_S152_m4_distalfemur/',
'data/Nuclei_and_CellsP40_S152_m4_proximaltibia/',
}; 


%embryo does not have height features so give wrong 
count=1;
for fi=[2,4,6]
    path=allpathEmbryo{fi};s=strsplit(path,'Nuclei_and_CellsE185_'); 
    dataFile1=strcat('../../EmbryoData/Femur_Tibia/MakeListNucleiLabelled/',s{2});
    d_femur{count}=load([dataFile1,'SaveVariablesForPlot_.mat']);
    count=count+1;
end

count=1;
for fi=[2,4,6]
    path=allpathPostnatal{fi};s=strsplit(path,'Nuclei_and_CellsP40_'); 
    dataFile2=strcat('../proximalTibia_DistalFemur/MakeListNucleiLabelled/',s{2});
    d_tibia{count}=load([dataFile2,'SaveVariablesForPlot_.mat']);
    count=count+1;
end


%legname={'0.75*c','c','1.25*c'};

legname={'E18.5 PT','P40 PT'};
mycolor={'ro-','bs-','g*-'};
facealpha=[0.3,0.2,0.15];


%d2{3}=load([dataFile,'SaveVariablesForPlot_larger.mat']);



% for i=1:size(d1.dataAvgMeanFeatures,1) % 19 properites 
%     for j=1:size(d1.dataAvgMeanFeatures,2) % xbin data 
%             data1=d1.dataAvgMeanFeatures{i,j};   % individual point 
%             data2=(d2.realization_AvgMeanFeatures{i}(j,:))';
%             %[size(data1), size(data2)]
%             [h,p]=ttest2(data1,data2);
%             pvaluetest{i}(j,1)=p;
%             
%     end
% end
% 




h=figure;
set(gcf, 'PaperSize', [5 3]); %7
set(gcf, 'PaperPosition', [0 0 5 3]);
for i=1:3
     df=d_femur{i};
     dt=d_tibia{i};
     df_data(:,i)=df.Ybin;
     dt_data(:,i)=dt.Ybin;
end

 v1=nanmean(df_data,2);  sd1=nanstd(df_data')';
 v2=nanmean(dt_data,2);  sd2=nanstd(dt_data')';

%  q(1)=plot(df.Xbin,v1,mycolor{1},'linewidth',1);  
%  hold on 
%  q(2)=plot(dt.Xbin,v2,mycolor{2},'linewidth',1);  
 

% index = find(~isnan(v1));  x=df.Xbin(index);
% min_y1=v1(index)-sd1(index); max_y1=v1(index)+sd1(index);  
% index=1:length(x);
% stacky2=(min_y1);stacky1=(max_y1);
% fillxx=x(index([1:end end:-1:1 1]));fillyy=[stacky1(index); stacky2(index(end:-1:1)); stacky1(1)];
% h=fill(fillxx,fillyy,mycolor{1}(1),'EdgeColor','none','facealpha',facealpha(1));
%                               
% index = find(~isnan(v2));  x=dt.Xbin(index);
% min_y1=v2(index)-sd2(index); max_y1=v2(index)+sd2(index);  
% index=1:length(x);
% stacky2=(min_y1);stacky1=(max_y1);
% fillxx=x(index([1:end end:-1:1 1]));fillyy=[stacky1(index); stacky2(index(end:-1:1)); stacky1(1)];
% h=fill(fillxx,fillyy,mycolor{2}(1),'EdgeColor','none','facealpha',facealpha(2));

    Feature.TitleName='';
    Feature.savefile=1;
    Feature.flag=0;
    Feature.Ylabel='<cluster density>';
    Feature.Legend=legname;
    Feature.Unit='';
    Feature.save_folder=dir2;
    Feature.SaveName='clusterDensity';
    statistical_test_function(df_data,dt_data,df.Xbin,Feature)

% legend(q,legname,'location','northeast')
% xlabel('Long axis of bone')
% ylabel('cluster density');
% saveas(h,[dir2,'clusterDensity.png'])
% close all 




for k=1:length(savename)
%     h=figure;
%     set(gcf, 'PaperSize', [5 3]); %7
%     set(gcf, 'PaperPosition', [0 0 5 3]);
    for i=1:3
           df=d_femur{i};
           dt=d_tibia{i};
           %q(i)=errorbar(d.Xbin(index),d.AvgMeanFeatures{k}(index),d.AvgStdFeatures{k}(index),mycolor{i},'linewidth',1);
           df_data(:,i)=df.AvgMeanFeatures{k};
           dt_data(:,i)=dt.AvgMeanFeatures{k};
           
    end
    
           
%            v1=nanmean(df_data,2);  sd1=nanstd(df_data')';
%            v2=nanmean(dt_data,2);  sd2=nanstd(dt_data')';
%            
%            q(1)=plot(df.Xbin,v1,mycolor{1},'linewidth',1);
%            hold on
%            q(2)=plot(dt.Xbin,v2,mycolor{2},'linewidth',1);
%            
%            
%            index = find(~isnan(v1));  x=df.Xbin(index);
%            min_y1=v1(index)-sd1(index); max_y1=v1(index)+sd1(index);  
%            index=1:length(x);
%            stacky2=(min_y1);stacky1=(max_y1);
%            fillxx=x(index([1:end end:-1:1 1]));fillyy=[stacky1(index); stacky2(index(end:-1:1)); stacky1(1)];
%            h=fill(fillxx,fillyy,mycolor{1}(1),'EdgeColor','none','facealpha',facealpha(1));
%            
%                     index = find(~isnan(v2));  x=dt.Xbin(index);
%                     min_y1=v2(index)-sd2(index); max_y1=v2(index)+sd2(index);  
%                     index=1:length(x);
%                     stacky2=(min_y1);stacky1=(max_y1);
%                     fillxx=x(index([1:end end:-1:1 1]));fillyy=[stacky1(index); stacky2(index(end:-1:1)); stacky1(1)];
%                     h=fill(fillxx,fillyy,mycolor{2}(1),'EdgeColor','none','facealpha',facealpha(2)); 
%     legend(q,legname,'location','best')  %'northeast'
%     
    
%    
%     if (k>=16)&(k<=18)
%         title('local','fontweight','normal')
%         %ylim([-0.05,1.06])
%         legend(q,legname,'location','south')  %'northeast'
%     end
%     
%      if (k>=19)&(k<=21)
%         title('global','fontweight','normal')
%         legend(q,legname,'location','south')  %'northeast'
%         %ylim([-0.69,1])
%      end
%       
    
     
%    
%     
%     xlim([0,1])
%     xlabel('Long axis of bone')
    flag=0;
    for ls=1:length(latexSymbol)
        if k==latexSymbol(ls)
            flag=1;
        end
    end
   
    
    
    Feature.TitleName=ylabelname{k};
    Feature.flag=flag;
    Feature.Ylabel=ylabelname{k};
    Feature.Legend=legname;
    Feature.Unit='';
    Feature.save_folder=dir2;
    Feature.SaveName=savename{k};
    if (k==26)|(k==27)|(k==28)|(k==37)|(k==35)|(k==29)|(k==30)|(k==31)|(k==32)|(k==22)|(k==23)|(k==24)|(k==25)|(k==1)|(k==3)|(k==38)|(k==13)
        Feature.savefile=0;
    else
        Feature.savefile=1;
        
        Feature.TitleName
        %dt_data
        statistical_test_function(df_data,dt_data,df.Xbin,Feature)
        
    end
    
    
            
%     saveas(h,[dir2,savename{k},'.png'])
%     close all 
end



% All graph properties 


















% 
% 
%  h=figure;
% XL=0.08;XR=0.02;XGap=0.05;Row=2;
% YT=0.06;YB=0.12;YGap=0.12;Col=4;
% Width=(1-XL-XR-((Col-1)*XGap))/Col;
% Height=(1-YT-YB-((Row-1)*YGap))/Row;
% YPos=1-YT-Height; 
% 
% set(gcf, 'PaperSize', [9 4]); %7
% set(gcf, 'PaperPosition', [0 0 9 4]);
% %mycolor={'b.-','g.-','c.-','r.-','k.-','m.-','b*-','r*-'};
% %mycolor={'r*-','r*-','r*-','r*-','r*-','r*-','r*-','r*-'};
% ylabelname={'3 tris','4 clique', '2 star',  '4 chordcycle', '4 tailed tris', '4 cycle', '3 star','4 path'};
% 
% for i=1:Row
%     XPos=XL;
%     for j=1:Col
%         chro=j+(i-1)*Col;
%             marray=[XPos,YPos,Width,Height];
%             subplot('Position',marray);
%             for k=1:length(dataFile)
%             
%                      load([dataFile{k},'SaveVariablesForPlot.mat'])
% 
%                      %p(k)=plot(Xbin,AvgMeanGraphlet{k},mycolor{k});
%                     plot(Xbin,AvgMeanGraphlet{chro}/max(AvgMeanGraphlet{chro}),mycolor{k},'linewidth',1);
%                     hold on 
%             %xlim([-0.05,1.05])
%             end
%             if i==2
%                 xlabel('Long axis of bone');
%             end
%             if j==1
%               ylabel('frequency/ maximum');
%             end
%            
%             title(ylabelname{chro},'fontweight','normal')
% 
% 	XPos=XPos+Width+XGap;
%     end
%     YPos=YPos-YGap-Height;
% end
% 
% %legend(p,ylabelname,'location','north');
% saveas(h,[dir2,'GraphletFrequency','.png'])
% close all





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



function output=dofitting(index,myinterval,celavg,pos1,pos2,pos3,mycolor,plotcolor)
               x=(myinterval(index))';
               y=celavg(index);
               output=errorAndPvalue(x,y); 
               q=plot(x,output.yfit,plotcolor,'linewidth',0.6);
%                text(x(1),pos1, ['Y = ',sprintf('%0.2f',output.p(1)),'X + ',sprintf('%0.2f',output.p(2))],'fontsize',5,'color',mycolor)
%                text(x(1),pos2, ['R^2 = ', sprintf('%0.2f',output.Rsq1)],'fontsize',5,'color',mycolor)
%                text(x(1),pos3, ['pvalue = ', sprintf('%0.2e',output.pvalue)],'fontsize',5,'color',mycolor)
               [output.Rsq1,output.Rsq2];
               output.plot=q;
end




function  output=errorAndPvalue(x,y)
               g=fittype(@(m,c,x) (m*x+c));              
               [p,s]=polyfit(x,y,1);
               [yfit,d] = polyval(p,x,s); 

               yresid = y-yfit;
               SSresid = sum(yresid.^2);
               SStotal = sum( (y-mean(y)).^2); % same thing   (length(y)-1) * var(y)
               nrsq = 1 - SSresid/SStotal;
               nrsq_adj = 1 - SSresid/SStotal * (length(y)-1)/(length(y)-length(p));
                 
               [f,stat]=fit(x,y,g,'startpoint',[-0.1,1]);
               
               % Two Sided 2*tcdf(2.29,99,'upper') 
               % One Sided   tcdf(2.29,99,'upper')
               %https://handbook-5-1.cochrane.org/chapter_7/7_7_3_2_obtaining_standard_deviations_from_standard_errors_and.htm
               %https://stattrek.com/regression/slope-test.aspx   
            
               SE = sqrt(SSresid/(length(y)-2)) / sqrt( sum((x-mean(x)).^2) ); % Standard Error 
               mypvalue=2*tcdf(abs(f.m/SE),stat.dfe, 'upper');
               [f.m/SE,stat.dfe,mypvalue];
               output.pvalue=mypvalue;
               output.yfit=yfit;
               output.Rsq1=nrsq;
               output.Rsq2=stat.rsquare;
               output.SE=SE;
               output.p=p;
end




function   boundary=statisticalTest(celavg,nucavg)
              significant=[];
              for ii=1:46
                  ind2=ii:ii+4;
                  [h,pP]=ttest2(celavg(ind2),nucavg(ind2));
                  if pP<0.05
                      significant=[significant,ii];
                  end
              end    
              
              boundary=[significant(1)];
              for ii=2:length(significant)
                  if significant(ii)-significant(ii-1)~=1
                      boundary=[boundary,significant(ii-1),significant(ii)];
                  end
              end
              boundary=[boundary,significant(end)+4];
end
