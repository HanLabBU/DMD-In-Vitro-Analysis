
% % in vivo data
%cd('\\engnas.bu.edu\research\eng_research_handata\EricLowet\DMD\invivoDMD\')
%cd('\\engnas.bu.edu\research\eng_research_handata\EricLowet\DMD\comp_wide_indi\')

close all;
clear all;

% Use the scripts in the DMD scripts folder
addpath('.');

% in vitro data
cd('~/handata_server/Pierre Fabris/DMD Project/All In Vitro Analysis/');

% Folder to save figures
%save_fig_path = '\\engnas.bu.edu\research\eng_research_handata\Pierre Fabris\DMD Project\Data Figures\';
save_fig_path = '~/hanlab_server/DMD Paper/In Vitro Plots/';

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
indiB=[];wideB=[]; indiSNR=[]; wideSNR=[]; indiAllB = []; wideAllB = [];
indi_SRate = []; wide_SRate = [];
indi_samp = []; wide_samp = [];
fov_label = {};
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
    
  if ignore_first == 1 && size(indifile.allresults.bleach, 1) > 1
    % Store the bleaching values except for the first trial
    
    % Ignore the second trial of 24 individual
    if ignore_24 == 1 && contains(indifile.allresults.fov_name, '24') == 1
        indi_temp = indifile.allresults.bleach(3:end, mROI(:, 1));
        wide_temp = widefile.allresults.bleach(3:end, mROI(:, 1));
    else
        indi_temp = indifile.allresults.bleach(2:end, mROI(:, 1));
        wide_temp = widefile.allresults.bleach(2:end, mROI(:, 1));
    end
    
  else
    
      if ignore_24 == 1 && contains(indifile.allresults.fov_name, '24') == 1
        % Ignore the second trial of the bleaching values
        indi_temp = indifile.allresults.bleach([1, 3:end], mROI(:, 1));
        wide_temp = widefile.allresults.bleach([1, 3:end], mROI(:, 1));
      else
        % Store all of the bleaching values
        indi_temp = indifile.allresults.bleach(:, mROI(:, 1));
        wide_temp = widefile.allresults.bleach(:, mROI(:, 1));
      end
  end

    % DEBUG
    if ~isempty(find(indi_temp < .5))
        disp(['Large Decay ' indifile.allresults.fov_name ' ' indifile.allresults.type]);
    end
    
    if ~isempty(find(wide_temp < .5))
        disp(['Large Decay ' widefile.allresults.fov_name ' ' widefile.allresults.type]);
    end
    
  
  indiAllB = horzcat_pad(indiAllB, indi_temp(:));
  wideAllB = horzcat_pad(wideAllB, wide_temp(:));
  
  % Store the FOV labels
  fov_label = cat(2, fov_label, [indifile.allresults.fov_name ' ' indifile.allresults.type]);
  fov_label = cat(2, fov_label, [widefile.allresults.fov_name ' ' widefile.allresults.type]);
  
  clear allIsnr
  for tr=1:size(  indifile.allresults.spike_snr,2)
      for ne=mROI(:,1)' %1:size(  indifile.allresults.spike_snr,1)
   allIsnr(ne,tr)= mean(indifile.allresults.spike_snr{ne,tr});
      end;end
    clear allWsnr
  for tr=1:size(  widefile.allresults.spike_snr,2)
      for ne= mROI(:,1)' %1:size(  widefile.allresults.spike_snr,1)
   allWsnr(ne,tr)= mean(widefile.allresults.spike_snr{ne,tr});
      end;end
    indiSNR= [indiSNR; nanmean(allIsnr,2) ];
    wideSNR= [wideSNR; nanmean(allWsnr,2) ];
  
    % Calculate the average event rate for each pattern
    indi_spikerate_neuron = [];
    wide_spikerate_neuron = [];
    
    for neuron=1:size(indifile.allresults.roaster, 1)
        total_spikes = sum(indifile.allresults.roaster(neuron, :));
        indi_spikerate_neuron = [indi_spikerate_neuron, total_spikes*sample_frequency./(size(indifile.allresults.roaster, 2))];
        total_spikes = sum(widefile.allresults.roaster(neuron, :));
        wide_spikerate_neuron = [wide_spikerate_neuron, total_spikes*sample_frequency./(size(widefile.allresults.roaster, 2))];
    end
    
    % Calculate the average spike event rate for the whole field of view
    indi_SRate = [indi_SRate, nanmean(indi_spikerate_neuron)];
    wide_SRate = [wide_SRate, nanmean(wide_spikerate_neuron)];
    
    % Store the spike amplitudes
    temp = [indifile.allresults.spike_amplitude{:}];
    indi_samp = horzcat_pad(indi_samp, temp(:));
    
    temp = [widefile.allresults.spike_amplitude{:}];
    wide_samp = horzcat_pad(wide_samp, temp(:));
    
end


%% PLOT the bleaching decay boxplots between DMD and Wide Field
figure('COlor','w'),plot(indiAllB(:),'r'); hold on,plot(wideAllB(:),'k')
legend indi wide
[h,p,ci,stats] = ttest(indiAllB(:),wideAllB(:))

%M=[ ((1-(indiB)).*-1)'  ,((1-(wideB)).*-1)'].*100;
figure('COlor','w','Renderer', 'painters')
boxplot([indiAllB(:), wideAllB(:)], {'Individual DMD', 'Wide Field'},  'notch','on',   'colors',[ 0.4 0.4 0.4], 'symbol','.k')
title_string = [ 'Boxplot of photodecay p= ' num2str(p)];
title(title_string);

saveas(gcf, [save_fig_path 'Photobleaching\Jpeg Format\' title_string '.jpg']);
saveas(gcf, [save_fig_path 'Photobleaching\EPS Format\' title_string '.eps'], 'epsc');
saveas(gcf, [save_fig_path 'Photobleaching\SVG Format\' title_string '.svg']);


% Violin plots of photobleaching
figure('Position', [300 300 800 750]);
violin(horzcat_pad(indiAllB(:), wideAllB(:)), 'xlabel', {'DMD', 'Wide Field'}, 'facecolor', [138/255 175/255 201/255]);
ylabel('Photobleach ratio');

title_string = [];
if ignore_first == 1
    title_string = ['Sumamry violin plots photobleaching ratios of individual DMD and wide field without first trials'];
else
    title_string = ['Summary violin plots photobleaching ratios of individual DMD and wide field'];
end
title(title_string);

saveas(gcf, [save_fig_path 'Photobleaching\Jpeg Format\' title_string '.jpg']);
saveas(gcf, [save_fig_path 'Photobleaching\EPS Format\' title_string '.eps'], 'epsc');

% Violin plot of photobleaching ratios by culture FOV
figure;
fov_bleach = [];
for i=1:size(indiAllB, 2)
    fov_bleach = horzcat_pad(fov_bleach, indiAllB(:, i));
    fov_bleach = horzcat_pad(fov_bleach, wideAllB(:, i));
end
violin(fov_bleach, 'xlabel', fov_label);
ax = gca;
xtickangle(ax, 45);
title_string = [];
if ignore_first == 1
    title_string = ['Photobleaching ratio DMD vs. Wide field by fov without first trials'];
else
    title_string = ['Photobleaching ratio DMD vs. Wide field by fov'];
end
title(title_string);
saveas(gcf, [save_fig_path 'Photobleaching\Jpeg Format\' title_string '.jpg']);
saveas(gcf, [save_fig_path 'Photobleaching\EPS Format\' title_string '.eps'], 'epsc');

% Plot of photobleaching ratios line graph between individual DMD and Wide
% Field
indi_pb_trial_ave = nanmean(indiAllB, 1);
wide_pb_trial_ave = nanmean(wideAllB, 1);
figure;
plot([indi_pb_trial_ave; wide_pb_trial_ave], '-', 'LineWidth', 3);
xticks([1, 2]);
xticklabels({'Individual DMD', 'Wide Field'});
xlim([0.5 2.5]);
title_string = [];
if ignore_first == 1
    title_string = ['Summary line plots photobleaching ratios of individual DMD and wide field without first trials'];
else
    title_string = ['Summary line plots photobleaching ratios of individual DMD and wide field'];
end
title(title_string);
saveas(gcf, [save_fig_path 'Photobleaching\Jpeg Format\' title_string '.jpg']);
saveas(gcf, [save_fig_path 'Photobleaching\EPS Format\' title_string '.eps'], 'epsc');

% Violin plot of the photo decay
figure('Position', [300 300 450 450]);
indi_decay = -100.*[repmat(1, length(indiAllB(:)), 1) - indiAllB(:)];
wide_decay = -100.*[ repmat(1, length(wideAllB(:)), 1) - wideAllB(:)];
violin(horzcat_pad(indi_decay, wide_decay), ...
    'xlabel', {'DMD', 'Wide Field'}, 'facecolor', [138/255 175/255 201/255]);

title_string = [];
if ignore_first == 1
    title_string = ['Summary violin plots photobleaching decay of individual DMD and wide field without first trials'];
else
    title_string = ['Summary violin plots photobleaching decay of individual DMD and wide field'];
end
title(title_string);
saveas(gcf, [save_fig_path 'Photobleaching\Jpeg Format\' title_string '.jpg']);
saveas(gcf, [save_fig_path 'Photobleaching\EPS Format\' title_string '.eps'], 'epsc');

% Plot the photobleaching decay bar graph 
% T test does not work here because the values are not correctly paired
figure('COlor','w','Position', [300 300 300 450]);
V1=nanmean(indi_decay);V1s=std(indi_decay)./sqrt(length(indi_decay));
V2=nanmean(wide_decay);V2s=std(wide_decay)./sqrt(length(wide_decay));
bar( [ 1 ], [V1],0.7,'FaceColor', [ 0.7 0.2 0.1]) , hold on,bar( [ 2 ], [V2],0.7,'FaceColor', [0.1 0.4 0.7])
set(gca,'Xtick', [ 1 2],'Xticklabel', {'DMD' ; 'Widefield'})
errorbar([ 1 2], [ V1 V2], [V1s V2s],'.k','Linewidth', 2)
axis tight;ylabel('Signal Reduction %')
title_string = [];
if ignore_first == 1
    title_string = ['Average signal decay without first trials'];
else
    title_string = ['Average signal decay'];
end
title(title_string);

saveas(gcf, [save_fig_path 'Photobleaching\Jpeg Format\' title_string '.jpg']);
saveas(gcf, [save_fig_path 'Photobleaching\EPS Format\' title_string '.eps'], 'epsc');


% Plot the boxplot SNRs
figure('COlor','w', 'Renderer', 'painters'),plot(indiSNR,'r'); hold on,plot(wideSNR,'k')
legend indi wide
[h,p,ci,stats] = ttest(indiSNR,wideSNR)

figure('COlor','w')
boxplot([indiSNR, wideSNR], {'Individual DMD', 'Wide Field'}, 'notch', 'on', 'colors',[ 0.4 0.4 0.4], 'symbol','k');
title_string = ['Boxplots of SNR p= ' num2str(p)];
title(title_string);

saveas(gcf, [save_fig_path 'SNR\Jpeg Format\' title_string '.jpg']);
saveas(gcf, [save_fig_path 'SNR\EPS Format\' title_string '.eps'], 'epsc');
saveas(gcf, [save_fig_path 'SNR\SVG Format\' title_string '.svg']);


%% Plot the number of resolvable event rate between individual and wide field

% Will be tricky because the event rate has to be over the total
% time of imaging and trials had different lengths
figure('Renderer', 'painters');
boxplot([indi_SRate', wide_SRate'], {'Individual DMD', 'Wide Field'}, 'notch', 'on', 'colors', [ 0.4 0.4 0.4], 'symbol','.k');
%hold on;
%plot([indi_SRate; wide_SRate], '--');
title_string = ['Resolved Spike Rate Individual DMD vs. Wide Field'];
ylabel('Spike Rate');
title(title_string);

saveas(gcf, [save_fig_path 'Event Rate\Jpeg Format\' title_string '.jpg']);
saveas(gcf, [save_fig_path 'Event Rate\EPS Format\' title_string '.eps'], 'epsc');
saveas(gcf, [save_fig_path 'Event Rate\SVG Format\' title_string '.svg']);

%% Plot the distribution of spike amplitudes
figure('Renderer', 'painters');
boxplot(horzcat_pad(indi_samp(:), wide_samp(:)), {'Individual DMD', 'Wide Field'}, 'notch', 'on', 'colors', [ 0.4 0.4 0.4], 'symbol','.k');
title_string = ['Spike amplitude Individual DMD vs. Wide Field'];
ylabel('delta F');
title(title_string);

saveas(gcf, [save_fig_path 'Spike Amplitude\Jpeg Format\' title_string '.jpg']);
saveas(gcf, [save_fig_path 'Spike Amplitude\SVG Format\' title_string '.svg']);