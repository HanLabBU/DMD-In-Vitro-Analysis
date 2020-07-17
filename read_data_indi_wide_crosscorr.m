
% In vivo path
%cd('\\engnas.bu.edu\research\eng_research_handata\EricLowet\DMD\invivoDMD\')
%cd('\\engnas.bu.edu\research\eng_research_handata\EricLowet\DMD\comp_wide_indi\')

close all;
clear all;

addpath('.');

% In vitro path
cd('\\engnas.bu.edu\research\eng_research_handata\Pierre Fabris\DMD Project\All In Vitro Analysis\');

% Final figure folder
save_fig_path = '\\engnas.bu.edu\research\eng_research_handata\Pierre Fabris\DMD Project\Data Figures\';

% Scripts for needed functions
addpath('\\ad\eng\research\eng_research_handata\EricLowet\DMD\main_analysis\');

ses=dir('*.mat');

%%% search for wide
clear findwide
for id=1:length(ses) 
 if   strfind(ses(id).name, 'Wide')>0
     findwide(id)=1;
 else
      findwide(id)=0;
 end
    
end
wideloc=find(findwide);

% Find corresponding individual field
indiloc = [];
for file=1:length(ses)
    if contains(ses(file).name, 'Individual') == 1
        load([ses(file).name]);
        indi_name = allresults.fov_name;

        for id=1:length(ses)
            if contains(ses(id).name, 'Wide') == 1 & contains(ses(id).name, indi_name) == 1           
                % Sanity check
                load([ses(id).name]);

                if strcmp(allresults.fov_name, indi_name) == 1
                    indiloc = [indiloc, file];
                end
            end
        end
    end
end

% Original code for finding the corresponding wide field mask
%wideloc=find(findwide);
%%%% search for fitting indi
%indiloc=[];
%for id2=1:length(wideloc)
%widename=ses(wideloc(id2)).name;
%    for id=1:length(ses) 
%        if id~=wideloc(id2)
%            indiname=   ses(id).name;
%            vx2= strfind(widename,'wide');;vx= strfind(indiname,'fov');
%            if length(indiname) > vx2(1)-1
%                if length(strfind(indiname(1:vx2(1)-1), widename(1:vx2(1)-1)))>0
%                    indiloc=[indiloc, id];
%                end
%            end
%        end
%    end
%end

 sc=0.1625*2; % micrometer per pixel 40x
%% seelct matching ROI
%indiB=[];wideB=[]; indiSNR=[]; wideSNR=[];
cindi=[];cwide=[]; rdist=[];
cindiS=[];cwideS=[];  allCx=[];allCm=[];
for id=1:length(wideloc)
    widefile=load(ses(wideloc(id)).name);
    indifile=load(ses(indiloc(id)).name);
    clear wrm % ROI centroid
    for id2= 1:length(widefile.allresults.roi)
        [ x y]=find(widefile.allresults.roi{id2});
        wrm(:,id2)= round(mean([x , y]));end
    clear irm
    for id2= 1:length(indifile.allresults.roi)
        [ x y]=find(indifile.allresults.roi{id2});
        irm(:,id2)= round(mean([x , y]));end  
    mROI=[];
    for id3=1:size(irm,2) % matching ROI
        cents=irm(:,id3);
        pxdiff=(sqrt(sum(bsxfun(@minus, wrm, cents).^2)));
    
      wloc=find(pxdiff<8);
      if ~isempty(wloc)
    mROI=[mROI; [ id3 wloc]];end
    end
    
    %%% Put in matrix%%%
 %subthresIndi= indifile.allresults.trace_ws;
  % subthresWide= widefile.allresults.trace_ws;  
    subthresIndi= indifile.allresults.orig_traceDN;
   subthresWide= widefile.allresults.orig_traceDN;

   % subthresIndi= subthresIndi-fastsmooth( subthresIndi,10,1,1);
    %  subthresWide= subthresWide-fastsmooth( subthresWide,10,1,1);
nsel=mROI(:,1);
clear allCmax  ROIdist
for ind1=1:length(nsel)
for ind2=1:length(nsel)
%  A1=  zscore(subthresIndi(nsel(ind1),:));A2=  zscore(subthresIndi(nsel(ind2),:));
   % [c,lags]= xcorr(A1,A2,10,'Coeff');
        A1=  (subthresIndi(nsel(ind1),:));A2=  (subthresIndi(nsel(ind2),:));
    [CC]=corrcoef(A1,A2, 'Rows', 'complete');
    allCmax(ind1,ind2)= CC(1,2);%max(abs(c));
   ROIdist(ind1,ind2)= sqrt(sum((irm(:,nsel(ind1)).*sc-irm(:,nsel(ind2)).*sc).^2));
end;end
    allCmax(allCmax>0.99)=NaN;
    
clear allCmaxW  
for ind1=1:length(nsel)
for ind2=1:length(nsel)
%  A1=  zscore(subthresWide(nsel(ind1),:));A2=  zscore(subthresWide(nsel(ind2),:));
    A1=  (subthresWide(nsel(ind1),:));A2=  (subthresWide(nsel(ind2),:));
 
   [CC]=corrcoef(A1,A2,'rows','complete');
    allCmaxW(ind1,ind2)= CC(1,2);%max(abs(c));
   ;
end;end
    allCmaxW(allCmaxW>0.99)=NaN;
      allCmaxW(allCmaxW==0)=NaN;
    
    cindi= [cindi; allCmax(:)];
cwide= [cwide; allCmaxW(:)];
rdist= [rdist; ROIdist(:)];

    %%
     subthresIndi= indifile.allresults.roaster;
   subthresWide= widefile.allresults.roaster;
rate_thres=5;
    clear allCmax  ROIdist
for ind1=1:length(nsel)
for ind2=1:length(nsel)
  A1=  zscore(subthresIndi(nsel(ind1),:));A2=  zscore(subthresIndi(nsel(ind2),:));
%  if ind1<ind2
 if  length(find(A1>0)) >rate_thres & length(find(A2>0)) >rate_thres
    [c,lags]= xcorr(fastsmooth(A1,10,1,1),fastsmooth(A2,10,1,1),100,'Coeff');
    c1=c(95:105);
    mean_corr = mean(c1);
    allCmax(ind1,ind2)= mean_corr;
    allCx= [allCx ; c];
    allCm=[allCm; mean_corr];
 %     allCx= [allCx ; c];
 %allCmax(ind2,ind1)= c1(n2);

 else;  allCmax(ind1,ind2)=NaN;end;
   ROIdist(ind1,ind2)= sqrt(sum((irm(:,nsel(ind1)).*sc-irm(:,nsel(ind2)).*sc).^2));
   % ROIdist(ind2,ind1)= ROIdist(ind1,ind2);
%  end
end;end
    allCmax(allCmax>0.99)=NaN;
      allCmax(allCmax==0)=NaN;
    
  clear allCmaxW  
for ind1=1:length(nsel)
for ind2=1:length(nsel)
  A1=  zscore(subthresWide(nsel(ind1),:));A2=  zscore(subthresWide(nsel(ind2),:));
    if ind1<ind2
   if  length(find(A1>0)) >rate_thres & length(find(A2>0)) >rate_thres
    [c,lags]= xcorr(fastsmooth(A1,10,1,1),fastsmooth(A2,10,1,1),100,'Coeff');
    c1=c(95:105) % This is the +/- around 0 from the 200 parameter, the offset
    mean_corr = mean(c1); % Use mean
    allCmaxW(ind1,ind2)= mean_corr;
     allCmaxW(ind2,ind1)= mean_corr;
   % allCmaxW(ind2,ind1)=  allCmax(ind1,ind2);
   else;  allCmaxW(ind1,ind2)=NaN;end
    end
end;end
    allCmaxW(allCmaxW>0.99)=NaN;  
         allCmaxW(allCmaxW==0)=NaN;
    
    
cindiS= [cindiS; allCmax(:)];
cwideS= [cwideS; allCmaxW(:)];


  
end
rdist2=rdist;
cindi=cindi(~isnan(cwide));
rdist=rdist(~isnan(cwide));
cwide=cwide(~isnan(cwide));

cindiS=cindiS(~isnan(cwideS))
%allCx=allCx(~isnan(cwideS),:)
rdist2=rdist2(~isnan(cwideS))
cwideS=cwideS(~isnan(cwideS))
rdist2=rdist2(~isnan(cindiS))
cwideS=cwideS(~isnan(cindiS))
cindiS=cindiS(~isnan(cindiS))
%allCx=allCx(~isnan(cindiS),:)

% PLOT subthreshold cross correlation
figure('Color','w')
plot(rdist,cindi,'.r','Markersize',20)
hold on, plot(rdist,cwide,'.k','Markersize',20)
fitResults1 = polyfit(rdist,cindi,1);
yplot1 = polyval(fitResults1,rdist);
plot(rdist,yplot1,'-r','Linewidth',2)
fitResults2 = polyfit(rdist,cwide,1);
yplot1 = polyval(fitResults2,rdist);
plot(rdist,yplot1,'-k','Linewidth',2)
xlabel('Distance')
ylabel('Corr')
legend({['indi slope=' num2str(fitResults1(1))], ['wide slope=' num2str(fitResults2(1))]});

title_string = ['Cross correlation subthreshold (Vm) indi vs. wide'];
title(title_string);
saveas(gcf, [save_fig_path 'Cross Correlation/Jpeg Format/' title_string '.jpg']);
saveas(gcf, [save_fig_path 'Cross Correlation/EPS Format/' title_string '.eps'], 'epsc');

[h,p,ci,stats] = ttest(cindi,cwide)
figure('COlor','w','Position', [ 300 300 200 200])
V1=nanmean(cindi);V1s=(std(cindi)./sqrt(length(cindi)));
V2=nanmean(cwide);V2s=(std(cindi)./sqrt(length(cwide)));
bar( [ 1 ], [V1],0.7,'FaceColor', [ 0.7 0.2 0.1]) , hold on,bar( [ 2 ], [V2],0.7,'FaceColor', [0.1 0.4 0.7])
set(gca,'Xtick', [ 1 2],'Xticklabel', {'DMD' ; 'Widefield'})
errorbar([ 1 2], [ V1 V2], [V1s V2s],'.k','Linewidth', 2)
axis tight;ylabel('corr max')
xlim([ 0.5 2.5]); %ylim([0  20])
title([ 'p= ' num2str(p)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spike to spike correlation


figure('Color','w')
plot(rdist2,cindiS,'.r','Markersize',20)
hold on,plot(rdist2,cwideS,'.k','Markersize',20)

fitResults1 = polyfit(rdist2,cindiS,1);
yplot1 = polyval(fitResults1,rdist2);
plot(rdist2,yplot1,'-r','Linewidth',2)
fitResults2 = polyfit(rdist2,cwideS,1);
yplot1 = polyval(fitResults2,rdist2);
plot(rdist2,yplot1,'-k','Linewidth',2)
xlabel('Distance');
ylabel('Corr');
legend({['indi slope=' fitResults1(1)], ['wide slope=' fitResults2(1)]});
title_string = ['Cross correlation spike-spike indi vs. wide'];
title(title_string);

saveas(gcf, [save_fig_path 'Cross Correlation/Jpeg Format/' title_string '.jpg']);
saveas(gcf, [save_fig_path 'Cross Correlation/EPS Format/' title_string '.eps'], 'epsc');

[h,p,ci,stats] = ttest(cindiS,cwideS)
figure('COlor','w','Position', [ 300 300 200 200])
V1=nanmean(cindiS);V1s=(std(cindiS)./sqrt(length(cindiS)));
V2=nanmean(cwideS);V2s=(std(cindiS)./sqrt(length(cwideS)));
bar( [ 1 ], [V1],0.7,'FaceColor', [ 0.7 0.2 0.1]) , hold on,bar( [ 2 ], [V2],0.7,'FaceColor', [0.1 0.4 0.7])
set(gca,'Xtick', [ 1 2],'Xticklabel', {'DMD' ; 'Widefield'})
errorbar([ 1 2], [ V1 V2], [V1s V2s],'.k','Linewidth', 2)
axis tight;ylabel('corr max')
xlim([ 0.5 2.5]); %ylim([0  20])
title([ 'p= ' num2str(p)])




% %% PLOT
% figure('COlor','w'),plot(indiB,'r'); hold on,plot(wideB,'k')
% legend indi wide
% [h,p,ci,stats] = ttest(indiB,wideB)
% figure('COlor','w','Position', [ 300 300 200 200])
% V1=((1-nanmean(indiB)).*-1).*100;V1s=(std(indiB)./sqrt(length(indiB))).*100;
% V2=((1-nanmean(wideB)).*-1).*100;V2s=(std(wideB)./sqrt(length(wideB))).*100;
% bar( [ 1 ], [V1],0.7,'FaceColor', [ 0.7 0.2 0.1]) , hold on,bar( [ 2 ], [V2],0.7,'FaceColor', [0.1 0.4 0.7])
% set(gca,'Xtick', [ 1 2],'Xticklabel', {'DMD' ; 'Widefield'})
% errorbar([ 1 2], [ V1 V2], [V1s V2s],'.k','Linewidth', 2)
% axis tight;ylabel('signal reduction %')
% xlim([ 0.5 2.5]); %ylim([0  20])
% title([ 'p= ' num2str(p)])



