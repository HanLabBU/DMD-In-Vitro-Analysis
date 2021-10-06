clear all

addpath(genpath('\\engnas.bu.edu\Research\eng_research_handata\EricLowet\'))
%%% Change path to the processed data (from data folder D)
path_to_data='\\engnas.bu.edu\research\eng_research_handata\EricLowet\DMD\comp_wide_indi\'
cd(path_to_data)

ses=dir('*.mat'); % FInd processed data

FS=500; % sampling rate
Noisefloor= 767.7 ; %% Camera noise floor measured experimentally 

%% search for wide
clear findwide
for id=1:length(ses)
    if   strfind(ses(id).name, 'wide')>0
        findwide(id)=1;
    else
        findwide(id)=0;
    end
    
end
wideloc=find(findwide);
%% search for fitting DMD
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



indiB=[];wideB=[]; indiSBR=[]; wideSBR=[];
 indirate=[]; widerate=[];
for id=1:length(wideloc)  %% main loop
    
    %% seelct matching ROI
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
    
    indiB= [indiB, nanmean(indifile.allresults.bleach(1:end,mROI(:,1)),1) ];
    wideB= [wideB, nanmean(widefile.allresults.bleach(1:end,mROI(:,1)),1)  ];
    MM=nanmean(indifile.allresults.bleach(1:end,mROI(:,1)),2);
    indiB_t(1:length(MM),id)= MM;
    MM=nanmean(widefile.allresults.bleach(1:end,mROI(:,1)),2);
    wideB_t(1:length(MM),id)= MM;
    
   
    
    clear allIsnr allIrate
    mm=0;
    for tr=1:size(  indifile.allresults.spike_snr,2)
        mm=mm+1;
        for ne=mROI(:,1)' 
            allIsbr(ne,mm)= nanmean(indifile.allresults.spike_snr{ne,tr});
            allIrate(ne,mm)= sum(indifile.allresults.roaster(ne,:))./(size(indifile.allresults.roaster,2)./FS);
        end;end
    
    clear allWsnr allWrate
    for tr=1:size(  widefile.allresults.spike_snr,2)
        for ne= mROI(:,1)' 
            allWsbr(ne,tr)= nanmean(widefile.allresults.spike_snr{ne,tr});
            allWrate(ne,tr)= sum(widefile.allresults.roaster(ne,:))./(size(widefile.allresults.roaster,2)./FS);
        end;end

    indiSBR= [indiSBR; nanmean(allIsbr,2) ];
    wideSBR= [wideSBR; nanmean(allWsbr,2) ];
    
    indirate= [indirate; nanmean(allIrate,2) ];
    widerate= [widerate; nanmean(allWrate,2) ];
    
end



%% PLOT PHOTOBLEACHING
M=[ (((indiB)).*-1)'  ,(((wideB)).*-1)'].*100
figure('COlor','w','Position',[ 300 300 250 200],'Renderer', 'painters')
boxplot( M   ,[ ones(length(indiB),1);ones(length(indiB),1).*2],  'notch','on',   'colors',[ 0.4 0.4 0.4], 'symbol','.k')
ylabel('Photobleaching %')
set(gca,'Xtick', [ 1 2],'Xticklabel', {'DMD' ; 'Widefield'})
title([ 'p= ' num2str(p)])

%% PLOT  SBR  
M=[ indiSBR  , wideSBR];
[h,p,ci,stats] = ttest(M(:,1),M(:,2))
figure('COlor','w','Position',[ 300 300 250 200],'Renderer', 'painters')
boxplot( M   ,[ ones(length(indiSBR),1);ones(length(wideSBR),1).*2],  'notch','on',   'colors',[ 0.4 0.4 0.4], 'symbol','.k')
set(gca,'Xtick', [ 1 2],'Xticklabel', {'DMD' ; 'Widefield'})
title([ 'p= ' num2str(p)])
ylabel('SBR')


%% PLOT FIRING RATE 
M=[indirate  , widerate];
[h,p,ci,stats] = ttest(M(:,1),M(:,2))
figure('COlor','w','Position',[ 300 300 250 200],'Renderer', 'painters')
boxplot( M   ,[ ones(length(indirate),1);ones(length(widerate),1).*2],  'notch','on',   'colors',[ 0.4 0.4 0.4], 'symbol','.k')
set(gca,'Xtick', [ 1 2],'Xticklabel', {'DMD' ; 'Widefield'})
title([ 'p= ' num2str(p)])
ylabel('Spike Hz')


