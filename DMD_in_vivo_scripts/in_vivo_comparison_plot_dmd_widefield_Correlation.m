clear all

addpath(genpath('\\engnas.bu.edu\Research\eng_research_handata\EricLowet\'))
%%% Change path to the processed data (from data folder D)
path_to_data='\\engnas.bu.edu\research\eng_research_handata\EricLowet\DMD\comp_wide_indi\'
cd(path_to_data)

ses=dir('*.mat'); % Find processed data

FS=500; % sampling rate

%%% search for widefield illumination recording data
clear findwide
for id=1:length(ses)
    if   strfind(ses(id).name, 'wide')>0
        findwide(id)=1;
    else
        findwide(id)=0;
    end
    
end
wideloc=find(findwide);
%%% search for DMD illumination recording data
indiloc=[];
for id2=1:length(wideloc)
    widename=ses(wideloc(id2)).name;
    for id=1:length(ses)
        if id~=wideloc(id2)
            indiname=   ses(id).name;
            vx2= strfind(widename,'wide');;vx= strfind(indiname,'fov');
            if length(indiname) > vx2(1)-1
                if length(strfind(indiname(1:vx2(1)-1), widename(1:vx2(1)-1)))>0
                    indiloc=[indiloc, id];
                end;end
        end
    end
end



sc=0.1625*2; % micrometer per pixel (40x objective)

%% seelct matching ROI
cindi=[];cwide=[]; rdist=[];
cindiS=[];cwideS=[];  allCx=[];allCm=[];

for id=1:length(wideloc) %% MAIN LOOP
    
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
    
    
    %% Correlation of subhtreshold somArchon traces (traces with removed spike segments)
    subthresIndi= indifile.allresults.trace_ws;
    subthresWide= widefile.allresults.trace_ws;
    
    nsel=mROI(:,1);
    if id==4
        nsel=mROI(1:4,1);  %session 612777_fov2_individual_allresults, ROI 4 and 5 the same
    end
    clear allCmax  ROIdist
    for ind1=1:length(nsel)
        for ind2=1:length(nsel)
            if ind1 < ind2
                A1=  (subthresIndi(nsel(ind1),:));A2=  (subthresIndi(nsel(ind2),:));
                [CC]=corrcoef(A1,A2,'Rows','complete');
                allCmax(ind1,ind2)= CC(1,2);%
                ROIdist(ind1,ind2)= sqrt(sum((irm(:,nsel(ind1)).*sc-irm(:,nsel(ind2)).*sc).^2));
            end;end;end
    allCmax(allCmax>0.99)=NaN;  %%  Autocorrelations are removed
    allCmax(allCmax==0)=NaN;    %%  empty fields in the matrix are excluded
    
    clear allCmaxW
    for ind1=1:length(nsel)
        for ind2=1:length(nsel)
            if ind1 < ind2
                A1=  (subthresWide(nsel(ind1),:));A2=  (subthresWide(nsel(ind2),:));
                [CC]=corrcoef(A1,A2,'Rows','complete');
                allCmaxW(ind1,ind2)= CC(1,2);%
                ;
            end;end;end
    allCmaxW(allCmax>0.99)=NaN;
    allCmaxW(allCmaxW==0)=NaN;
    
    cindi= [cindi; allCmax(:)];
    cwide= [cwide; allCmaxW(:)];
    rdist= [rdist; ROIdist(:)];
    
    %% Correlation of spike trains
    subthresIndi= indifile.allresults.roaster;
    subthresWide= widefile.allresults.roaster;
    rate_thres=0;
    sm= 10; % spike train smoothinger kernel , in 2ms bins
    clear allCmax  ROIdist
    for ind1=1:length(nsel)
        for ind2=1:length(nsel)
            if ind1 < ind2
                A1=  zscore(subthresIndi(nsel(ind1),:));A2=  zscore(subthresIndi(nsel(ind2),:));
                if  length(find(A1>0)) >=rate_thres & length(find(A2>0)) >=rate_thres
                    [c,lags]= xcorr(fastsmooth(A1,sm,1,1),fastsmooth(A2,sm,1,1),200,'Coeff');
                    c1=c(196:206); % +- 8ms
                    [n1 n2]= max(abs(c1));
                    allCmax(ind1,ind2)= mean(c1);
                    allCx= [allCx ; c];
                    allCm=[allCm;c1(n2)];
                    
                else;  allCmax(ind1,ind2)=NaN;end;
                ROIdist(ind1,ind2)= sqrt(sum((irm(:,nsel(ind1)).*sc-irm(:,nsel(ind2)).*sc).^2));
            end;end;end
    allCmax(allCmax>0.99)=NaN;
    allCmax(allCmax==0)=NaN;
    
    clear allCmaxW   allCmaxW
    for ind1=1:length(nsel)
        for ind2=1:length(nsel)
            if ind1<ind2
                A1=  zscore(subthresWide(nsel(ind1),:));A2=  zscore(subthresWide(nsel(ind2),:));
                
                if  length(find(A1>0)) >=rate_thres & length(find(A2>0)) >=rate_thres
                    [c,lags]= xcorr(fastsmooth(A1,sm,1,1),fastsmooth(A2,sm,1,1),200,'Coeff');
                    c1=c(196:206);% +- 8ms
                    [n1 n2]= max(abs(c1));
                    allCmaxW(ind1,ind2)= mean(c1);
                else;  allCmaxW(ind1,ind2)=NaN;end
            end
        end;end
    allCmaxW(allCmaxW>0.99)=NaN;
    allCmaxW(allCmaxW==0)=NaN;
    
    
    cindiS= [cindiS; allCmax(:)];
    cwideS= [cwideS; allCmaxW(:)];
    
    
    
end
rdist2=rdist;
% Remove any NaNs
cindi=cindi(~isnan(cwide));
rdist=rdist(~isnan(cwide));
cwide=cwide(~isnan(cwide));
cindiS=cindiS(~isnan(cwideS))
rdist2=rdist2(~isnan(cwideS))
cwideS=cwideS(~isnan(cwideS))
rdist2=rdist2(~isnan(cindiS))
cwideS=cwideS(~isnan(cindiS))
cindiS=cindiS(~isnan(cindiS))


%% PLOT subthreshold Vm correlation
[h,p,ci,stats] = ttest(cindi,cwide)
M=[cindi;cwide];
figure('COlor','w','Position',[ 300 300 250 200],'Renderer', 'painters')
boxplot( M   ,[ ones(length(cindi),1);ones(length(cwide),1).*2],  'notch','on',   'colors',[ 0.4 0.4 0.4], 'symbol','*k')
ylabel('Corr')
set(gca,'Xtick', [ 1 2],'Xticklabel', {'DMD' ; 'Widefield'})
title([ 'p= ' num2str(p)])

%% PLOT spike train correlation
[h,p,ci,stats] = ttest(cindiS,cwideS)
M=[cindiS;cwideS];
fig=figure('COlor','w','Position',[ 300 300 250 200],'Renderer', 'painters');%set(gca, 'FontName', 'Arial')
boxplot( M   ,[ ones(length(cindiS),1);ones(length(cwideS),1).*2],  'notch','on',   'colors',[ 0.4 0.4 0.4], 'symbol','*k')
ylabel('Corr');%set(fig,'DefaultAxesFontName','Arial')
set(gca,'Xtick', [ 1 2],'Xticklabel', {'DMD' ; 'Widefield'})
title([ 'p= ' num2str(p)]);%set(fig, 'FontName', 'Arial')
ylim([-0.05 0.2])


