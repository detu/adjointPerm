% input
 workDir        = [ROOTDIR 'TEST1\'];
 time           = 10950;
 iter           = 1;

% production data file
 summaryFile    = [workDir 'summary_' num2str(time) '_iter' num2str(iter) '.mat'];
 load(summaryFile,'wellData','wellDataE','timeSteps','settings');
 
 ntimes = length(timeSteps);
 nwells = length(wellData);
 nr = 4; nc = 3; fig =0;
 
% Well rate time series 
 ip = 0; fig = fig + 1; figure(fig);
 for iw = 1 : nwells
     
     rate(:,iw) = wellData(iw).flx;
     wct(:,iw) = wellData(iw).ffl;     
     if exist('wellDataE','var')
       rateE(:,:,iw) = wellDataE(iw).flx;
       wcte(:,:,iw) = wellDataE(iw).ffl;       
       ne = size(rateE,2);
     else
       ne = 0;
     end
     
     figure(1)
     ip = ip + 1;
     subplot(nr,nc,ip)
     % ensemble members
     if ne > 0
       mprodE = zeros(length(timeSteps),1);  
       for j=1:ne
         prodE = squeeze(rateE(:,j,iw));
         prodE = prodE * day; % m^3/s -> m^3/day          
         h = plot(timeSteps,prodE,'-k'); cl=[0.5 0.5 0.5]; hold on
         set(h,'Linestyle','-','Linewidth',1,'Color',cl,'Marker','none','MarkerSize',2,'MarkerFaceColor',cl,'MarkerEdgeColor',cl);
         mprodE=mprodE+prodE./ne;
       end  
       % ensemble mean 
       h = plot(timeSteps,mprodE,'-k'); cl='b';
       set(h,'Linestyle','-','Linewidth',1,'Color',cl,'Marker','none','MarkerSize',2,'MarkerFaceColor',cl,'MarkerEdgeColor',cl);
     end
     % truth
     prod = squeeze(rate(:,iw));
     prod = prod * day; % m^3/s -> m^3/day
     h = plot(timeSteps,prod,'-k'); cl='r';
     set(h,'Linestyle','-','Linewidth',1,'Color',cl,'Marker','none','MarkerSize',2,'MarkerFaceColor',cl,'MarkerEdgeColor',cl);
         
     set(gca,'XLim',[0 timeSteps(end)]); %,'YLim',[ymin ymax])
     set(gca,'Fontsize',11,'Linewidth',1); 
     title(['Q (m^3/day) well ' num2str(iw)],'Fontsize',11);
           
 end
 
 % Water cut time series
 ip = 0; fig = fig + 1; figure(fig);
 for iw = 1 : nwells
     
     wct(:,iw) = wellData(iw).ffl;     
     if exist('wellDataE','var')
       wctE(:,:,iw) = wellDataE(iw).ffl;       
       ne = size(wctE,2);
     else
       ne = 0;
     end
     
     figure(2)
     ip = ip + 1;
     subplot(nr,nc,ip)
     % ensemble members
     if ne > 0
       mprodE = zeros(length(timeSteps),1);  
       for j=1:ne
         prodE = squeeze(wctE(:,j,iw));        
         h = plot(timeSteps,prodE,'-k'); cl=[0.5 0.5 0.5]; hold on
         set(h,'Linestyle','-','Linewidth',1,'Color',cl,'Marker','none','MarkerSize',2,'MarkerFaceColor',cl,'MarkerEdgeColor',cl);
         mprodE=mprodE+prodE./ne;
       end  
       % ensemble mean 
       h = plot(timeSteps,mprodE,'-k'); cl='b';
       set(h,'Linestyle','-','Linewidth',1,'Color',cl,'Marker','none','MarkerSize',2,'MarkerFaceColor',cl,'MarkerEdgeColor',cl);
     end
     % truth
     prod = squeeze(wct(:,iw));
     h = plot(timeSteps,prod,'-k'); cl='r';
     set(h,'Linestyle','-','Linewidth',1,'Color',cl,'Marker','none','MarkerSize',2,'MarkerFaceColor',cl,'MarkerEdgeColor',cl);
         
     set(gca,'XLim',[0 timeSteps(end)],'YLim',[0 1]);
     set(gca,'Fontsize',11,'Linewidth',1); 
     title(['WCT well ' num2str(iw)],'Fontsize',11);
     
 end

 % Pressure time series
 ip = 0; fig = fig + 1; figure(fig);
 for iw = 1 : nwells
     
     bhp(:,iw) = wellData(iw).bhp;
     if exist('wellDataE','var')
       bhpE(:,:,iw) = wellDataE(iw).bhp;       
       ne = size(bhpE,2);
     else
       ne = 0;
     end
     
     figure(3)
     ip = ip + 1;
     subplot(nr,nc,ip)
     % ensemble members
     if ne > 0
       mprfE = zeros(length(timeSteps),1);         
       for j=1:ne
         prfE = squeeze(bhpE(:,j,iw));
         prfE = prfE / barsa;         
         h=plot(timeSteps,prfE,'-k'); cl=[0.5 0.5 0.5]; hold on
         set(h,'Linestyle','-','Linewidth',1,'Color',cl,'Marker','none','MarkerSize',2,'MarkerFaceColor',cl,'MarkerEdgeColor',cl);
         mprfE = mprfE + prfE./ne;
       end
       % ensemble mean
       h = plot(timeSteps,mprfE,'-k'); cl='b';
       set(h,'Linestyle','-','Linewidth',1,'Color',cl,'Marker','none','MarkerSize',2,'MarkerFaceColor',cl,'MarkerEdgeColor',cl);
     end
     % truth
     prf = squeeze(bhp(:,iw));
     prf = prf / barsa;
     h = plot(timeSteps,prf,'-k'); cl='r';
     set(h,'Linestyle','-','Linewidth',1,'Color',cl,'Marker','none','MarkerSize',2,'MarkerFaceColor',cl,'MarkerEdgeColor',cl);
     
     set(gca,'XLim',[0 timeSteps(end)]); %,'YLim',[ymin ymax])
     set(gca,'Fontsize',11,'Linewidth',1);
     title(['P (bar) well ' num2str(iw)],'Fontsize',11);
     
 end
 
%% printing 
%  set(gcf,'PaperUnits','normalized')
%  set(gcf,'PaperType','a4letter')
%  set(gcf,'PaperPosition',[.01 .01 .99 .99])
%  set(gcf,'PaperOrientation','Landscape')
%  figname=fullfile(pwd,['plotProduction.pdf']);
%  print(gcf,'-dpdf',figname);
