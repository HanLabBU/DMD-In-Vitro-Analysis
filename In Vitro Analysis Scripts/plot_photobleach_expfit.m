close all;
clc;

%% Author notes for using script
% Will be using the trial{1-3}.traces to do the spike detection

t_start = tic();

% Save figure directory
save_fig_path = '~/handata_server/Pierre Fabris/DMD Project/Data Figures/';

addpath('.');
addpath('~/handata_server/EricLowet/DMD/main_analysis');

% Data directory
invitro_data_path = '~/handata_server/Pierre Fabris/DMD Project/All In Vitro Analysis/';
cd(invitro_data_path);

ses=dir(['*.mat']);


% Find all of the 
% Search for wide
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
        indifile = load([ses(file).name]);
        indi_name = indifile.allresults.fov_name;

        for id=1:length(ses)
            if contains(ses(id).name, 'Wide') == 1 & contains(ses(id).name, indi_name) == 1           
                % Sanity check
                widefile = load([ses(id).name]);
                if strcmp(widefile.allresults.fov_name, indi_name) == 1 & strcmp(widefile.allresults.type, 'wide') == 1
                    indiloc = [indiloc, file];
                end
            end
        end
    end
end

% Sanity check that matching pairs are found
if length(indiloc) ~= length(wideloc)
	disp('Not all of the data found are matched pairs');
	pause;
end

% Store the SNRs and total number of spikes with the corresponding threshold, so {<threshold>, <SNRs>}
indi_SNR = [];
wide_SNR = [];

% Loop through all of the data files
for id=1:length(indiloc) 
    
    % Printing id
    id
   
	indifile = load(ses(indiloc(id)).name);
	widefile = load(ses(wideloc(id)).name);
	indiresults = indifile.allresults;
	wideresults = widefile.allresults;

    % Find all of the matching ROIs for corresponding 
    mROI = []; % [<corresponding neurons>, <indi roi id, matched wide roi id>]
    
    indiCentroids = [];
    for roi_id = 1:length(indifile.allresults.roi)
        [x y] = find(indifile.allresults.roi{roi_id});
        indiCentroids = [indiCentroids; round(mean([x, y], 1))];    
    end
    
    wideCentroids = [];
    for roi_id = 1:length(widefile.allresults.roi)
        [x y] = find(widefile.allresults.roi{roi_id});
        wideCentroids = [wideCentroids; round(mean([x, y], 1))];    
    end
    
    for roi_id = 1:size(indiCentroids)
        indiCentroid = indiCentroids(roi_id, :);
        pxdist = sqrt(sum( bsxfun(@minus, wideCentroids, indiCentroid).^2, 2) );
        wideROILoc = find(pxdist < 8);
        
        if ~isempty(wideROILoc)
            mROI = [mROI; roi_id, wideROILoc];
        end
    end

    % Can use the traces_photo_fit function from Eric's folder
    % Loop through each ROI from the individual
    
    % [<fit parameter>, <neuron>]
    indi_photofit = [];

    for trace_id=mROI(:, 1)'
        %Store the average of photobleach values across trials for each neuron
        trial_photo_vals = [];
        
        for trial_id=1:length(indifile.allresults.trial)
            cur_trace = indifile.allresults.trial{trial_id}.traces(:, trace_id);
            filtered_trace = medfilt1(cur_trace, 51);
            [fitval, fitgood, fitgood2] = traces_photo_fit(filtered_trace - 767.7); % Floor subtracted
        end

    end
    
end


%% Plot the SNR
figure('Renderer', 'painters');
indi_err = nanstd(indi_SNR, 0, 2);
wide_err = nanstd(wide_SNR, 0, 2);
errorbar(up_thres', nanmedian(indi_SNR, 2), indi_err, '--r');
hold on;
errorbar(up_thres',nanmedian(wide_SNR, 2), wide_err, '--b');
ylabel('Spike SBR');
xlabel('Spike SBR Threshold Value');
title_string = 'Detection Threshold vs. SNR';
legend({'Targeted', 'Widefield'});
title(title_string);

% Print the 4.5 median values
%indi_med = nanmedian(indi_SNR(find(up_thres == 4.5), :), 2)
%wide_med = nanmedian(wide_SNR(find(up_thres == 4.5), :), 2)


%DEBUG
total_indi_SNRS = length(indi_SNR)
total_wide_SNRS = length(wide_SNR)
total_ROIs

% Perform paired t-test of all the SNRs between the conditions for each threshold
% Save everything in a matrix to be saved to a csv file
stats_table = ["Threshold", "Median Targeted", "Median Widefield", "P-Value", "Confidence interval", "Degrees of Freedom", "T-statistic", "sd (indi - wide)"];
for i=1:length(up_thres)
    disp(['Paired-sample t-test ' num2str(up_thres(i))]);
    [h, p, ci, stats] = ttest(indi_SNR(i, :), wide_SNR(i, :))
    stats_table = [stats_table; ...
                   up_thres(i), nanmedian(indi_SNR(i, :)), nanmedian(wide_SNR(i, :)), p, string(num2str(ci)), stats.df, stats.tstat, stats.sd];
end

% Save the stats table
stats_table
writematrix(stats_table, [save_fig_path 'Vary Thres/threshold vs. SNR statistics table.csv']);

% Print the final time
t_final = toc(t_start)
