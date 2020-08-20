%% Set the root path of all the data
close all;
clear all;

% In vitro data network
root_path = '\\ad\eng\research\eng_research_handata\Pierre Fabris\DMD Project\v2 script Trace Extraction\In vitro\';

% In vitro local folder
%root_path = 'D:\Work\Graduate School\Han Lab\Projects\Digital_Micromirror_Device\Temporary In Vitro Culture Data\';

% One folder that shows both the individual mask and wide field
%root_path = '\\ad\eng\research\eng_research_handata\Yangyang Wang\DMD invivo data processing\July 4 processing by folder\D\602088 40x obj individual vs large fov 1\';

% Folder to save all analysis data
save_all_path = '\\ad\eng\research\eng_research_handata\Pierre Fabris\DMD Project\All In Vitro Analysis\';

% Exclude trials with obvious motion artefacts
exclude_motart = 1;
artefact_trials = {{'Culture 5\Wide', 2}, {'Culture 24\Wide', 3}, {'Culture 23\Wide', 3}}; % {Session, <trial no>}

% Store all of the directories that have trace data
% Depending on how the data is organized, this script will perform analyses
% for all trials in a given folder

% Recursively find all FOV directories from the root path
FOV_dirs = {};
search_dirs = {root_path};
while ~isempty(search_dirs)
    % Pop the first directory of the list
    cur_dir = search_dirs{1};
    search_dirs;
    
    % Search through directory and find more directories
    listings = dir(cur_dir);
    listings = listings(~ismember({listings.name},{'.','..'}));
    dirs = {listings.name};
    isdirs = [listings.isdir];

    dir_exists = 0;
    for i=1:length(dirs)
        % Add directories to search_dirs
        if isdirs(i)
            search_dirs{end + 1} = [cur_dir dirs{i} '\'];
            dir_exists = 1;
        end
    end
    
    % Current directory is an FOV folder if no subdirectories were found
    if ~dir_exists
        FOV_dirs{end + 1} = cur_dir;
    end
    
    % Remove first element which was the cur_dir
    search_dirs(1) = [];
end

%% Save the ROIs, the original traces, the name of the FOV folder, and
% the mask type
for i = 1:length(FOV_dirs)
    fov = FOV_dirs{i}
    listings = dir(fov);
    files = {listings.name};
    
    % Check if all results file exists for this FOV in the saved folder
    allresults_file = [strrep(erase(fov, root_path), '\', '_') 'allresults']
    allresults.fov_name = strrep(erase(fov, {root_path, 'Individual Mask', 'Wide Field'}), '\', '');
    
    % Read in all trials and save original traces
    trace_files = files(contains(files, 'traces'));
    allresults.trial = {};
    num_skipped = 0;
    
    for j=1:length(trace_files)
        
        % Check if this trial is listed to have motion artefacts
        skip_trial = 0;
        if exclude_motart
            for k=1:length(artefact_trials)
                if contains(fov, artefact_trials{k}{1}) == 1 & j == artefact_trials{k}{2}
                    skip_trial = 1;
                end
            end
        end
        
        if skip_trial == 1
            % DEBUG
            [fov ' trial # ' num2str(j)]
            
            disp('Skipped');
            num_skipped = num_skipped + 1;
            continue
            
        end
        
        load([fov trace_files{j}]);
        allresults.trial{j - num_skipped}.traces = traces;
        
	end
    
    % Save all the ROIs
    roi_file = files(contains(files, 'ROIs'));
    load([fov roi_file{1}]);
    allresults.roi = ROIs;
    
    % Save the type of mask
    if contains(fov, 'Individual Mask')
        allresults.type = 'indi';
    elseif contains(fov, 'Wide Field')
        allresults.type = 'wide';
    end
    allresults
    save([save_all_path allresults_file], 'allresults');
end

%% Calculate the cross correlations and distance between neurons
% TODO fix to iterate through allresults
fov_results = dir([save_all_path '*.mat']);
length({fov_results.name})
for i = 1:length({fov_results.name})
    
    % Check if allresults file exists for this FOV in the saved folder
    allresults_file = fov_results(i).name
    load([save_all_path allresults_file]);
    
    % Compile all of the traces in the trial
    max_trial = 0;
    for j=1:length(allresults.trial)
        max_trial = max(max_trial, size(allresults.trial{j}.traces, 1));
    end
    
    trial_traces = NaN(max_trial, size(allresults.trial{1}.traces, 2));
    for j=1:length(allresults.trial)
        trial_traces(1:size(allresults.trial{j}.traces, 1), :, j) = allresults.trial{j}.traces;
    end
    
    cross_corr_result = cross_correlation_distance(trial_traces, allresults.roi);
    
    allresults.roi_centroids = cross_corr_result.centroids;
    allresults.neuron_dist = cross_corr_result.neuron_dist;
    
    save([save_all_path allresults_file], 'allresults');
end

%% Detect spikes and calculate SNR
% At the moment, different tweaks are going to be made to the SNR algorithm
fov_results = dir([save_all_path '*.mat']);
length({fov_results.name})
for i = 1:length({fov_results.name})
    allresults_file = {fov_results.name};
    allresults_file = allresults_file{i}
    load([save_all_path allresults_file]);
    
    for j=1:length(allresults.trial)
        traces = allresults.trial{j}.traces;
        
        % Perform the spike_detection_SNR
        result = spike_detect_SNR_v3b(traces); % v3 was the original script used for the previous
        
        % Store each trial's spike and SNR data columnwise
        if j == 1
            allresults.orig_trace = result.orig_trace;
            allresults.denoise_trace = result.denoise_trace;
            allresults.trace_ws = result.trace_ws;
            allresults.orig_traceDN = result.orig_traceDN;
            allresults.roaster = result.roaster;
            allresults.spike_snr = result.spike_snr;
            allresults.spike_amplitude = result.spike_amplitude;
            allresults.spike_idx = result.spike_idx;
            allresults.trace_noise = result.trace_noise;        
        else
            % Save all of the SNR results data into the allresults structure
            allresults.orig_trace = horzcat(allresults.orig_trace, result.orig_trace);
            allresults.denoise_trace = horzcat(allresults.denoise_trace, result.denoise_trace);
            allresults.trace_ws = horzcat(allresults.trace_ws, result.trace_ws);
            allresults.orig_traceDN = horzcat( allresults.orig_traceDN, result.orig_traceDN);
            allresults.roaster = horzcat( allresults.roaster, result.roaster);
            allresults.spike_snr = horzcat(allresults.spike_snr, result.spike_snr);
            allresults.spike_amplitude = horzcat(allresults.spike_amplitude, result.spike_amplitude);
            allresults.spike_idx = horzcat( allresults.spike_idx, result.spike_idx);
            allresults.trace_noise = horzcat(allresults.trace_noise, result.trace_noise);        
        end
    end
    
    % Save all the SNR results 
    save([save_all_path allresults_file], 'allresults');
end

%% Calculate the photobleaching ratios for each FOV folder
% Store all of the photobleaching ratios
all_pb_ratios = [];
pb_plot_labels = {};
fov_results = dir([save_all_path '*.mat']);
for i = 1:length({fov_results.name})
    load([save_all_path fov_results(i).name]);
	fov_results(i).name
    
    % Store all of the ratios in a single folder
    folder_pb_ratios = [];
    
    % Trial mean/std of photobleaching ratio
    trial_pb_averages_std = [];
    
    % Load each trial's traces
    for j=1:length(allresults.trial)
            traces = allresults.trial{j}.traces;
            trial_pb_ratios = photobleach_absolute_values(traces, 300);
            folder_pb_ratios = [folder_pb_ratios; trial_pb_ratios];
            
            trial_pb_averages_std = [trial_pb_averages_std, [nanmean(trial_pb_ratios); std(trial_pb_ratios)]];
    end
    
    % Average and std of pb ratio for entire folder
    average_folder_pb = [nanmean(folder_pb_ratios); std(folder_pb_ratios)];
    
    % Save the photobleaching ratios into the all results variable
    allresults.bleach = folder_pb_ratios;
    
    % Add folder ratios to all pb ratio matrix
    all_pb_ratios = horzcat_pad(all_pb_ratios, folder_pb_ratios(:));
    
    [allresults.fov_name ' ' allresults.type]
    pb_plot_labels = cat(2, pb_plot_labels, [allresults.fov_name ' ' allresults.type]);
    
    % Update/Save all results
    save([save_all_path fov_results(i).name], 'allresults');
end


% Plot all of the folder photobleachings
figure;
boxplot(all_pb_ratios, pb_plot_labels);
ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
xtickangle(ax, 45);
ylabel('(initial intensity)/(final intensity)');


