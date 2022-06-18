clear all 


allpath={
'data/Nuclei_and_CellsP40_S151_m2_distalfemur/',
'data/Nuclei_and_CellsP40_S151_m2_proximaltibia/',
'data/Nuclei_and_CellsP40_S152_m3_distalfemur/',
'data/Nuclei_and_CellsP40_S152_m3_proximaltibia/',
'data/Nuclei_and_CellsP40_S152_m4_distalfemur/',
'data/Nuclei_and_CellsP40_S152_m4_proximaltibia/',
}; 


% Vdf=[];
% for gi=[1,3,5]
%     path=allpath{gi};
%     s=strsplit(path,'Nuclei_and_CellsP40_');
%      out2=strcat('MakeListNucleiLabelled/',s{2});    
%      a=load([out2,'AllFeaturesSave_.mat']);
%      Vdf=[Vdf;a.normalized_volume];
% end
% 
% Vpt=[];
% for gi=[2,4,6]
%     path=allpath{gi};
%     s=strsplit(path,'Nuclei_and_CellsP40_');
%      out2=strcat('MakeListNucleiLabelled/',s{2});    
%      a=load([out2,'AllFeaturesSave_.mat']);
%      Vpt=[Vpt;a.normalized_volume];
% end
% 
% Vdf=Vdf(find(Vdf<10));
% Vpt=Vpt(find(Vpt<10));
% 
% 
% subplot(1,2,1)
% hist(Vdf,30);
% subplot(1,2,2)
% hist(Vpt,30)





for gi=1:length(allpath)
		path=allpath{gi};
		disp(path)
        s=strsplit(path,'Nuclei_and_CellsP40_');
        out1=strcat('MakeListNuclei_unlabelled/',s{2});    
        out2=strcat('MakeListNucleiLabelled/',s{2});    
        out3=strcat('MakeListClustersMask/',s{2});  
        a1=load([out1,'centroid_and_surface_nuclei.mat']);
        a2=load([out2,'centroid_and_surface_nuclei.mat']);
        a3=load([out3,'centroid_and_surface_nuclei.mat']);
        
        c1=a1.centroid - mean(a1.centroid);
        c2=a2.centroid - mean(a2.centroid);
        c3=a3.centroid - mean(a3.centroid);
        ankit{gi}(1,:)=[min(c1(:,3)) max(c1(:,3))];
        ankit{gi}(2,:)=[min(c2(:,3)) max(c2(:,3))];
        ankit{gi}(3,:)=[min(c3(:,3)) max(c3(:,3))];

        h=figure;
        plot3(a1.centroid(:,1),a1.centroid(:,2),a1.centroid(:,3),'b.','markersize',3);
        hold on 
        plot3(a2.centroid(:,1),a2.centroid(:,2),a2.centroid(:,3),'r.','markersize',3);
        plot3(a3.centroid(:,1),a3.centroid(:,2),a3.centroid(:,3),'k.','markersize',3);
        xlabel('x');ylabel('y');zlabel('z')

        sname=s{2}(1:strlength(s{2})-1);
        saveas(h,['matlabFigures/', sname,'.fig'])
        saveas(h,['matlabFigures/',sname,'.png']);
        close all 
        

end
        
