
% In vivo path
%cd('\\engnas.bu.edu\research\eng_research_handata\EricLowet\DMD\invivoDMD\')
%cd('\\engnas.bu.edu\research\eng_research_handata\EricLowet\DMD\comp_wide_indi\')

close all;
clear all;

addpath('.');

% In vitro path
cd('~/handata_server/Pierre Fabris/DMD Project/All In Vitro Analysis/');
%cd('D:\DMD Analysis Temp Data\');


% Scripts for needed functions
addpath('~/handata_server/EricLowet/DMD/main_analysis/');

ses=dir('*.mat');

%%% search for wide
clear findwide
for id=1:length(ses) 
    if strfind(ses(id).name, 'Wide')>0
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

num_match_rois = [];
% Loop through each session
for id=1:length(wideloc)
    widefile=load(ses(wideloc(id)).name);
    indifile=load(ses(indiloc(id)).name);
    

    %% seelct matching ROI
    clear wrm % ROI centroid
    for id2= 1:length(widefile.allresults.roi)
        [ x y]=find(widefile.allresults.roi{id2});
        wrm(:,id2)= round(mean([x , y]));
    end
    clear irm    
    for id2= 1:length(indifile.allresults.roi)
        [ x y]=find(indifile.allresults.roi{id2});
        irm(:,id2)= round(mean([x , y]));
    end  
    mROI=[];
    
    for id3=1:size(irm,2) % matching ROI
        cents=irm(:,id3);
        pxdiff=(sqrt(sum(bsxfun(@minus, wrm, cents).^2)));
        wloc=find(pxdiff<8);
    
        if ~isempty(wloc)
            mROI=[mROI; [ id3 wloc]];
        end
    end
    %% End of ROI matching
    
    num_match_rois = [num_match_rois, length(mROI(:, 1)) ];    

end

num_match_rois

mean = nanmean(num_match_rois)
stdev = nanstd(num_match_rois)
sem = stdev/(sqrt(length(num_match_rois)))
