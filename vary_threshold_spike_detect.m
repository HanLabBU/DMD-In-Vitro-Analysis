clear all;
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
			% [snr vals, neuron]
			neuron_snrs = [];
			
			% Iterating through each neuron's SNRs
			for i=1:length(temp_var)
				x = [temp_var{i, 1}'];
                		if isempty(x), x = [NaN]; end;
				neuron_snrs = horzcat_pad(neuron_snrs, x');
			end
			
			thres_SNR = horzcat_pad(thres_SNR, neuron_snrs');

        		% Total spikes
			trial_num_spikes = [trial_num_spikes, sum(sum(result.roaster))];
		end
			
		% [neuron, threshold]
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
			
			% Iterating through each neuron's SNRs
			for i=1:length(temp_var)
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
	
	%TODO I think I should just list all SNRs into
	%indi_SNR = cat(3, indi_SNR, indi_fov_SNR);
	%indi_total_spikes = cat(3, );
	
	% Store everything as one matrix (no individual FOVs are discernible)
	indi_SNR = horzcat_pad(indi_SNR, indi_fov_SNR');
%	indi_total_spikes = [indi_total_spikes, indi_fov_num_spikes];

	wide_SNR = horzcat_pad(wide_SNR, wide_fov_SNR');
%	wide_total_spikes = [wide_total_spikes, wide_fov_num_spikes];
	
end


%% Plot the SNR
figure;
indi_err = nanstd(indi_SNR, 0, 2);
wide_err = nanstd(wide_SNR, 0, 2);
errorbar(up_thres', nanmedian(indi_SNR, 2), indi_err, '--r');
hold on;
errorbar(up_thres',nanmedian(wide_SNR, 2), wide_err, '--b');
ylabel('SNRs');
xlabel('Threshold Value');
title_string = 'Detection Threshold vs. SNR';
legend({'Indi', 'Wide'});
title(title_string);

% DEBUG
% Print the 4.5 median values
indi_med = nanmedian(indi_SNR(find(up_thres == 4.5), :), 2)
wide_med = nanmedian(wide_SNR(find(up_thres == 4.5), :), 2)

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
