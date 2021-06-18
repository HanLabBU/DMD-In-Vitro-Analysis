
% % in vivo data
%cd('\\engnas.bu.edu\research\eng_research_handata\EricLowet\DMD\invivoDMD\')
%cd('\\engnas.bu.edu\research\eng_research_handata\EricLowet\DMD\comp_wide_indi\')
clc;
close all;
clear all;
 
% Use the scripts in the DMD scripts folder
addpath('.');

% in vitro data
cd('~/handata_server/Pierre Fabris/DMD Project/All In Vitro Analysis/');

% Folder to save figures
%save_fig_path = '\\engnas.bu.edu\research\eng_research_handata\Pierre Fabris\DMD Project\Data Fi226 - 67 (NaNs) = gures\';
save_fig_path = '~/handata_server/Pierre Fabris/DMD Project/Data Figures/';

% Ignore the first trial across each FOV for photobleaching estimation
ignore_first = 1;
% Ignore the 2nd trial of 24
ignore_24 = 1;

% Original Sampling Frequency
sample_frequency = 500; % Hz

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
                if strcmp(allresults.fov_name, indi_name) == 1 & strcmp(allresults.type, 'wide') == 1
                    indiloc = [indiloc, file];
                    ['Individual ' ses(file).name]
                    ['Wide Field ' ses(id).name]
                end
            end
        end
    end
end

%% seelct matching ROI
indiSNR=[]; wideSNR=[]; indiAllB = []; wideAllB = [];
indi_SRate = []; wide_SRate = [];
indi_samp = []; wide_samp = [];
fov_label = {};
indi_no_spike = [];
wide_no_spike = [];
indi_Bdiff = [];
wide_Bdiff = [];
indi_Bdiff_norm = [];

for id=1:length(wideloc)
    indi_temp = [];
    wide_temp = [];
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

    % -- Bleaching data starts here

  if ignore_first == 1 & size(indifile.allresults.bleach, 1) > 1
    % Store the bleaching values except for the first trial
    
    % Ignore the second trial of 24 individual
    if ignore_24 == 1 & contains(indifile.allresults.fov_name, '24') == 1
        indi_temp = indifile.allresults.bleach(3:end, mROI(:, 1));
        wide_temp = widefile.allresults.bleach(3:end, mROI(:, 1));
    

    else
        indi_temp = indifile.allresults.bleach(2:end, mROI(:, 1));
        wide_temp = widefile.allresults.bleach(2:end, mROI(:, 1));
        

        % Find the decays that are weird
        [low_r, low_c] = find(indi_temp < 0.90);
        
        % Get all of the more "normal" decay
        [high_r, high_c] = find(indi_temp >= 0.90);

        % Loop through the lower decay
        for iter=1:length(low_r)
            
            indi_raw_trace = indifile.allresults.trial{low_r(iter)}.traces(:, low_c(iter));
            wide_raw_trace = widefile.allresults.trial{low_r(iter)}.traces(:, low_c(iter));
            
            figure;
            subplot(2, 2, 1);
            plot(indi_raw_trace);
            title('Raw trace');
            subplot(2, 2, 2);
            plot(medfilt1(indi_raw_trace, 51));
            title('Median filtered');
            subplot(2, 2, 3);
            plot(wide_raw_trace);
            subplot(2, 2, 4);
            plot(medfilt1(wide_raw_trace, 51));
            sgtitle([ses(indiloc(id)).name ' Trial ' num2str(low_r(iter)) ' Neuron ' num2str(low_c(iter))]);
            
            % Store trace difference from filtered part
            filtered = medfilt1(indi_raw_trace, 51);
            low_int = nanmean(filtered(end-300:end));
            high_int = nanmean(filtered(1:300));
            indi_Bdiff = [indi_Bdiff, high_int - low_int];
        
            filtered = medfilt1(wide_raw_trace, 51);
            low_int = nanmean(filtered(end-300:end));
            high_int = nanmean(filtered(1:300));
            wide_Bdiff = [wide_Bdiff, high_int - low_int];
        end
        
        for iter=1:length(high_r)
            
            indi_raw_trace = indifile.allresults.trial{high_r(iter)}.traces(:, high_c(iter));
            
            %figure;
            %subplot(2, 2, 1);
            %plot(indi_raw_trace);
            %title('Raw trace');
            %subplot(2, 2, 2);
            %plot(medfilt1(indi_raw_trace, 51));
            %title('Median filtered');
            %subplot(2, 2, 3);
            %plot(wide_raw_trace);
            %subplot(2, 2, 4);
            %plot(medfilt1(wide_raw_trace, 51));
            %sgtitle([ses(indiloc(id)).name ' Trial ' num2str(low_r(iter)) ' Neuron ' num2str(low_c(iter))]);
            
            % Store trace difference from filtered part
            filtered = medfilt1(indi_raw_trace, 51);
            low_int = nanmean(filtered(end-300:end));
            high_int = nanmean(filtered(1:300));
            indi_Bdiff_norm = [indi_Bdiff_norm, high_int - low_int];
        
        end
    end
    
  else
    
      if ignore_24 == 1 & contains(indifile.allresults.fov_name, '24') == 1
        % Ignore the second trial of the bleaching values
        indi_temp = indifile.allresults.bleach([1, 3:end], mROI(:, 1));
        wide_temp = widefile.allresults.bleach([1, 3:end], mROI(:, 1));
      else
        % Store all of the bleaching values
        indi_temp = indifile.allresults.bleach(:, mROI(:, 1));
        wide_temp = widefile.allresults.bleach(:, mROI(:, 1));
      end
  end
  
    %DEBUG figure out what is going on with the Fluorescence decay plot
    if size(indi_temp < 0.90) > 1
        disp(['Large bleach trial ' ses(indiloc(id)).name]);
    end

  indi_temp = (1 - indi_temp).*-100;   %./indi_tr_lengths;
  wide_temp = (1 - wide_temp).*-100;    %./wide_tr_lengths;
  
  indiAllB = horzcat_pad(indiAllB, nanmean(indi_temp,1));
  wideAllB = horzcat_pad(wideAllB, nanmean(wide_temp,1));
 
end


%% PLOT the bleaching decay boxplots between DMD and Wide Field
% figure('COlor','w'),plot(indiAllB(:),'r'); hold on,plot(wideAllB(:),'k')
% legend indi wide
% title('Photodecay line plot');

%M=[ ((1-(indiB)).*-1)'  ,((1-(wideB)).*-1)'].*100;
figure();
M=[ indi_Bdiff', wide_Bdiff'];
boxplot(M, {'Individual DMD', 'Wide Field'},  'notch','on',   'colors',[ 0.4 0.4 0.4], 'symbol','.k');
title('Raw difference between initial and end part of trace');

figure();
M=[ horzcat_pad(indi_Bdiff', indi_Bdiff_norm')];
boxplot(M, {'DMD wack decay', 'DMD normal decay'},  'notch','on',   'colors',[ 0.4 0.4 0.4], 'symbol','.k');
title('Raw difference');
