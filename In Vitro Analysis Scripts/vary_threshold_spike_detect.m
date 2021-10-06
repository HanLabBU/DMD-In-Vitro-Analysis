close all;
clc;

%% Author notes for using script
% Will be using the trial{1-3}.traces to do the spike detection

t_start = tic();

% Save figure directory
save_fig_path = '~/handata_server/Pierre Fabris/DMD Project/Data Figures/';

addpath('.');

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
indi_total_spikes = [];
wide_total_spikes = [];


% Range of threshold values to test
up_thres = [3:0.5:6];

% DEBUG
% This will be used to see if the number of indi_SNR elements is the same as the # of SNR values in the plots
total_ROIs = 0;


% Loop through all of the data files
for id=1:length(indiloc) 
    
    % Printing id
    id
   
	indifile = load(ses(indiloc(id)).name);
	widefile = load(ses(wideloc(id)).name);
	indiresults = indifile.allresults;
	wideresults = widefile.allresults;


	indi_fov_SNR = [];
	indi_fov_num_spikes = [];
	wide_fov_SNR = [];
	wide_fov_num_spikes = [];

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
    
    % Keep track of total ROIs used
    total_ROIs = total_ROIs + size(mROI, 1);

    %DEBUG
    mROI
    

	% Calculate values for the individual DMD
	for thres = up_thres
		trial_num_spikes = [];
		
		% Store each threshold SNR as [neuron, trial snrs]
		thres_SNR = [];
		for trial_id = 2:length(indiresults.trial)
 		   	
			% SNR
			result = spike_detect_SNR_v3b(indiresults.trial{trial_id}.traces, thres);	
			temp_var = result.spike_snr;
			
			% Extract and store each neuron's SNR into a matrix
			% Row, Col array representation: [<snr vals>, <neuron>]
			neuron_snrs = [];
			
			% Iterating through each neuron's SNRs with matched individual ROI
			for i=mROI(:, 1)'
				x = [temp_var{i, 1}'];
                		if isempty(x), x = [NaN]; end;
				neuron_snrs = horzcat_pad(neuron_snrs, x');
			end
			
            % Row, Col array representation: [<SNR at ROI>, <for trial #>]
			thres_SNR = horzcat_pad(thres_SNR, neuron_snrs');

        		% Total spikes
			trial_num_spikes = [trial_num_spikes, sum(sum(result.roaster))];
		end
			
		% [<neuron>, <SNR at given threshold>
        	indi_fov_SNR = horzcat_pad(indi_fov_SNR, nanmean(thres_SNR, 2));
		indi_fov_num_spikes = horzcat_pad(indi_fov_num_spikes, trial_num_spikes');
	end

	% Calculate values for wide field 
	for thres = up_thres
		trial_num_spikes = [];
		
		% Store each threshold SNR as [neuron, trial snrs]
		thres_SNR = [];
		for trial_id = 2:length(wideresults.trial)
 		   	
			% SNR
			result = spike_detect_SNR_v3b(wideresults.trial{trial_id}.traces, thres);	
			temp_var = result.spike_snr;
			
			% Extract and store each neuron's SNR into a matrix
			% [snr vals, neuron]
			neuron_snrs = [];
			
			% Iterating through each neuron's SNRs at the matched widefield ROI
            for i=mROI(:, 2)'
				x = [temp_var{i, 1}'];
                		if isempty(x), x = [NaN]; end;
				neuron_snrs = horzcat_pad(neuron_snrs, x');	
			end
			
			thres_SNR = horzcat_pad(thres_SNR, neuron_snrs');

        		% Total spikes
			trial_num_spikes = [trial_num_spikes, sum(sum(result.roaster))];
		end
			
		% [neuron, threshold]
        	wide_fov_SNR = horzcat_pad(wide_fov_SNR, nanmean(thres_SNR, 2));
		wide_fov_num_spikes = horzcat_pad(wide_fov_num_spikes, trial_num_spikes');  
	end
	
	%TODO I think I should just list all individual spike SNRs into the matrix. Ok, I gave up because it was a little complicated
	%indi_SNR = cat(3, indi_SNR, indi_fov_SNR);
	%indi_total_spikes = cat(3, );
	
	% Store everything as one matrix (no individual FOVs are discernible)
	indi_SNR = horzcat_pad(indi_SNR, indi_fov_SNR');
%	indi_total_spikes = [indi_total_spikes, indi_fov_num_spikes];

	wide_SNR = horzcat_pad(wide_SNR, wide_fov_SNR');
%	wide_total_spikes = [wide_total_spikes, wide_fov_num_spikes];
	
end

%% Plot the SNR as a line plot
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

% Create X-axis labels
xaxis_labels = repmat(["Targeted", "Widefield"], 1, length(up_thres));

% Interleave all of the SNRs for boxplotting
all_snrs = [];
for i=1:length(up_thres)
    
    % Append the threshold value to the label
    xaxis_labels(i*2 - 1) = xaxis_labels(i*2 - 1) + " " + num2str(up_thres(i));
    xaxis_labels(i*2) = xaxis_labels(i*2) + " " + num2str(up_thres(i));
    
    % Matrix is becoming [<SNRs from each FOV>, <threshold>]
    all_snrs = horzcat_pad(all_snrs, indi_SNR(i, :)');
    all_snrs = horzcat_pad(all_snrs, wide_SNR(i, :)');

end

%% Plot the SNRs as boxplots
figure('Renderer', 'painters');
boxplot(all_snrs, 'labels', xaxis_labels, 'notch', 'on', 'colors', [ 0.4 0.4 0.4], 'symbol', '.k');
title_string = 'Detection Threshold vs. SNR Boxplots';
title(title_string);

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

% % Plot the total number of spikes
% figure;
% plot(up_thres, sum(indi_total_spikes, 2), '--r');
% hold on;
% plot(up_thres, sum(wide_total_spikes, 2), '--b');
% ylabel('Number of spikes');
% xlabel('Threshold Value');
% legend({'Indi', 'Wide'});
% title_string = 'Detection Threshold vs. Number of Spikes';
% title(title_string);
% 

% Print the final time
t_final = toc(t_start)
