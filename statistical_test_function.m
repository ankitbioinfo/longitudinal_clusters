function []=statistical_test_function(cel,nuc,myinterval,Feature)

titlename=Feature.TitleName;
titlename='';
ylabelfeature=Feature.Ylabel;
FigureLegend=Feature.Legend;
unit=Feature.Unit;

          
%celallcolor={'r-','b-','g-','m-','k-'};
mycolor={'r','b','g','m','c'};
lt={'.','*','x','+','s'};

fcelallcolor={'ro','ro','ro','ro','ro','ro'};
fnucallcolor={'bs','bs','bs','b^','b^','b^'};

tname={'wt', 'mut'};     



h1=figure();
XL=0.07;XR=0.03;XGap=0.07;Row=1;
YT=0.06;YB=0.11;YGap=0.15;Col=2;
Width=(1-XL-XR-((Col-1)*XGap))/Col;
Height=(1-YT-YB-((Row-1)*YGap))/Row;
YPos=1-YT-Height; 

set(gcf, 'PaperSize', [10 5]);
set(gcf, 'PaperPosition', [0 0 10 5]);

% 1 for std deviaton test and 2 for mean test 
testwithstd=1;


for i=1:Row
    XPos=XL;
    for j=1:Col
        chro=j+(i-1)*Col;
        if chro<=6
            marray=[XPos,YPos,Width,Height];
            subplot('Position',marray);
            testwithstd=chro;
           


        celavg=nanmean(cel,2);
        nucavg=nanmean(nuc,2);
        
        celstd=(nanstd(cel'))';
        nucstd=(nanstd(nuc'))';
       
        p(1)=plot(myinterval,celavg,strcat(fcelallcolor{1},'-'),'linewidth',1);
        hold on 
        p(2)=plot(myinterval,nucavg,strcat(fnucallcolor{1},'-'),'linewidth',1);
        
        
        index = find(~isnan(celavg));  
        min_y1=celavg(index)-celstd(index); max_y1=celavg(index)+celstd(index);       
        x=myinterval;stacky2=(min_y1)';stacky1=(max_y1)';
        fillxx=x(index([1:end end:-1:1 1]));fillyy=[stacky1(index) stacky2(index(end:-1:1)) stacky1(1)];
        h=fill(fillxx,fillyy,'r','EdgeColor','none','facealpha',0.3);
        
        index = find(~isnan(nucavg));  
        min_y2=nucavg(index)-nucstd(index); max_y2=nucavg(index)+nucstd(index);  
        stacky2=(min_y2)';stacky1=(max_y2)';
        fillxx=x(index([1:end end:-1:1 1]));fillyy=[stacky1(index) stacky2(index(end:-1:1)) stacky1(1)];
        
     
        h=fill(fillxx,fillyy,'b','EdgeColor','none','facealpha',0.2);

        globalmin=min([min(min_y1),min(min_y2)]);
        globalmax=max([max(max_y1),max(max_y2)]);
      
 

            set(gca,'fontsize',11);            
            gap=[globalmin,globalmax];
            axis([0,1, globalmin, 1.2*globalmax]);

            
            if chro==2
            title(['mean test: ',titlename],'fontweight','normal','fontsize',12)
            else
            title(['std dev test: ',titlename],'fontweight','normal','fontsize',12)
            end

            
            n=length(myinterval);

            if testwithstd==1
              boundary=statisticalTest(celstd,nucstd);
              howmanybound=0;
              for ii=1:2:length(boundary)
                  %disp('bone');
                  ind2=boundary(ii):boundary(ii+1);
                  if length(ind2)>3
                          [h,pR]=ttest2(celstd(ind2),nucstd(ind2)); %stat= mes(celavg(ind2),nucavg(ind2),'hedgesg' );
                       if pR<0.05 
                          if mod(howmanybound,2)==0
                              value=1.1*globalmax;
                          else
                              value=1.1*globalmax-0.09*diff(gap);
                          end
                          if pR<0.001
                              text(ind2(1)/n,value,strcat('p=',sprintf('%0.1e',pR)));
                          else
                              text(ind2(1)/n,value,strcat('p=',sprintf('%0.3f',pR)));
                          end
                          howmanybound=howmanybound+1;
              %text(ind2(1)/50,value-0.09*diff(gap),strcat('g=',sprintf('%0.3f',abs(stat.hedgesg))));
              plot([ind2(1)/n,ind2(1)/n],[globalmin,2*globalmax],'k:','linewidth',0.5);
              plot([ind2(end)/n,ind2(end)/n],[globalmin,2*globalmax],'k:','linewidth',0.5);
              %[ind2(1)/50,ind2(end)/50]
              %rectangle('Position',[ind2(1)/50,0, (ind2(end)-ind2(1))/50, 11 ],'FaceColor',[0 .5 .5],'EdgeColor','none','LineWidth',0.1,'facealpha',0.02);
              
              x=[ind2(1)/n,ind2(end)/n];index=[1,2];stacky1=[globalmin,globalmin];stacky2=[2*globalmax,2*globalmax];
              fillxx=x(index([1:end end:-1:1 1]));fillyy=[stacky1(index) stacky2(index(end:-1:1)) stacky1(1)];
              h=fill(fillxx,fillyy,'k','EdgeColor','none','facealpha',0.05);
                       end
                  end
              end
            end
              
              
             if testwithstd==2  
              boundary=statisticalTest(celavg,nucavg);
              howmanybound=0;
              for ii=1:2:length(boundary)
                  %disp('bone');
                  ind2=boundary(ii):boundary(ii+1);
                  if length(ind2)>3
                          [h,pR]=ttest2(celavg(ind2),nucavg(ind2)); %stat= mes(celavg(ind2),nucavg(ind2),'hedgesg' );
                       if pR<0.05 
                          if mod(howmanybound,2)==0
                              value=1.1*globalmax;
                          else
                              value=1.1*globalmax-0.09*diff(gap);
                          end
                          if pR<0.001
                              text(ind2(1)/n,value,strcat('p=',sprintf('%0.1e',pR)));
                          else
                              text(ind2(1)/n,value,strcat('p=',sprintf('%0.3f',pR)));
                          end
                          howmanybound=howmanybound+1;
              %text(ind2(1)/50,value-0.09*diff(gap),strcat('g=',sprintf('%0.3f',abs(stat.hedgesg))));
              plot([ind2(1)/n,ind2(1)/n],[globalmin,2*globalmax],'k:','linewidth',0.1);
              plot([ind2(end)/n,ind2(end)/n],[globalmin,2*globalmax],'k:','linewidth',0.1);
              %[ind2(1)/50,ind2(end)/50]
              %rectangle('Position',[ind2(1)/50,0, (ind2(end)-ind2(1))/50, 11 ],'FaceColor',[0 .5 .5],'EdgeColor','none','LineWidth',0.1,'facealpha',0.02);
              
              x=[ind2(1)/n,ind2(end)/n];index=[1,2];stacky1=[globalmin,globalmin];stacky2=[2*globalmax,2*globalmax];
              fillxx=x(index([1:end end:-1:1 1]));fillyy=[stacky1(index) stacky2(index(end:-1:1)) stacky1(1)];
              h=fill(fillxx,fillyy,'k','EdgeColor','none','facealpha',0.05);
                       end
                  end
              end
             end
              
            
          
            if j==1
                %ylabel([ylabelfeature,' ',unit]);
                if Feature.flag==1
                    ylabel([ylabelfeature,' ',unit],'Interpreter','Latex');
                else
                    ylabel([ylabelfeature,' ',unit]);
                end
            end
            
            
          
            xlabel('Long axis of bone')
            legend(p,FigureLegend,'location','north','fontsize',7);
            legend 'boxoff'; 
            %legend(q,legendarray,'location','northwest','fontsize',7);

        
    
      end
	XPos=XPos+Width+XGap;
    end
    YPos=YPos-YGap-Height;
end


            if ~exist([Feature.save_folder],'dir')
               mkdir([Feature.save_folder]);
            end
            
            if Feature.savefile==1
            saveas(h1,[Feature.save_folder, Feature.SaveName]);
            saveas(h1,[Feature.save_folder, Feature.SaveName,'.png']);
            end
            close all 


end



function   boundary=statisticalTest(celavg,nucavg)
              n=size(celavg,1);
              significant=[];
              for ii=1:n-4
                  ind2=ii:ii+4;
                  [h,pP]=ttest2(celavg(ind2),nucavg(ind2));
                  if pP<0.05
                      significant=[significant,ii];
                  end
              end    
              
              if length(significant)>1
                  boundary=[significant(1)];
              else
                  boundary=[];
              end
              for ii=2:length(significant)
                  if significant(ii)-significant(ii-1)~=1
                      boundary=[boundary,significant(ii-1),significant(ii)];
                  end
              end
              if length(significant)>1
                boundary=[boundary,significant(end)+4];
              end
end