%% Set the root path of all the data
close all;
clear all;

% In vitro data network
root_path = '~/handata_server/Pierre Fabris/DMD Project/v2 script Trace Extraction/In vitro/';

% Folder to save all analysis datak
save_all_path = '~/handata_server/Pierre Fabris/DMD Project/All In Vitro Analysis/';

% Exclude trials with obvious motion artefacts
exclude_motart = 1;
artefact_trials = {{'Culture 5/Wide', 2}, {'Culture 24/Wide', 3}, {'Culture 23/Wide', 3}, {'Culture 29/Wide', 1}, {'Culture 29/Indi', 1}, {'Culture 24/Indi', 2} }; % {Session, <trial no>}

% Exclude weird looking traces
weird_traces = {{'1', 'Individual Mask', 3, [1, 2, 4]},  {'1', 'Wide Field', 1, [2, 6, 12]},  ...
{'1', 'Wide Field', 2, [1, 2, 4, 15]}, ...
{'21', 'Individual Mask', 2, 2}, {'21', 'Wide Field', 1, [2, 3]}, {'21', 'Wide Field', 2, 2}, ...
{'23', 'Individual Mask', 1, 2}, {'23', 'Individual Mask', 3, 5}, {'23', 'Wide Field', 2, 6}, ...
{'24', 'Individual Mask', 1, 13}, {'24', 'Individual Mask', 3, [10, 11]}, {'24', 'Individual Mask', 4, [10, 11]}, {'24', 'Wide Field', 1, 18}, {'24', 'Wide Field', 3, 10}, ...
{'25', 'Individual Mask', 1, [4, 6]}, {'25', 'Individual Mask', 3, 7}, {'25', 'Wide Field', 1, 6}, {'25', 'Wide Field', 2, 1}, {'25', 'Wide Field', 4, [2, 3, 4]}, ...
{'26', 'Individual Mask', 3, 4}, {'26', 'Individual Mask', 4, [2, 4, 18]}, ...
{'26', 'Wide Field', 1, [6, 7]}, {'26', 'Wide Field', 3, [2, 4, 18]},  {'26', 'Wide Field', 4, 18}, ...
{'3', 'Wide Field', 1, [4, 9, 11]}, ...
{'4', 'Individual Mask', 2, [3, 4]}, {'4', 'Individual Mask', 3, [3, 8]}, {'4', 'Wide Field', 2, 4}, ...
}; % {Session, condition, <trial no>, <neuron no.>};

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
            search_dirs{end + 1} = [cur_dir dirs{i} '/'];
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
    allresults_file = [strrep(erase(fov, root_path), '/', '_') 'allresults']
    allresults.fov_name = strrep(erase(fov, {root_path, 'Individual Mask', 'Wide Field'}), '/', '');
    
    % Read in all trials and save original traces
    trace_files = files(contains(files, 'traces'));
    allresults.trial = {};
    num_skipped = 0;
    
    % Loop through each trial file
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
        
        % Load this trial's traces
        load([fov trace_files{j}]);
        
        % TODO Remove individual traces here
        for k=1:length(weird_traces)
            if strcmp(erase(allresults.fov_name, 'Culture '), weird_traces{k}{1}) == 1 & ...
                    contains(fov, weird_traces{k}{2}) & j == weird_traces{k}{3}
                
                % Remove the trace
                traces(:, weird_traces{k}{4}) = [NaN];
                
                % DEBUG
                disp(['Removed ' allresults.fov_name ' ' fov ' Trial ' num2str(j) ' traces']);
                
            end
        end
        
        % TODO maybe the best option would be to have an empty array for trials that were skipped
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
    
    % Check to make sure trials are stored
    if length(allresults.trial) == 0
        trial_traces = [];
    else
        trial_traces = NaN(max_trial, size(allresults.trial{1}.traces, 2));
    end
    
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
        result = spike_detect_SNR_v3b(traces, 4.5); % v3 was the original script used for the previous
        
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
            allresults.orig_trace_untrended = result.orig_trace_untrended;
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
            allresults.orig_trace_untrended = horzcat(allresults.orig_trace_untrended, result.orig_trace_untrended);
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
	fov_results(i).name;
    
    % Store all of the ratios in a single folder
    folder_pb_ratios = [];
    
    % Load each trial's traces and calculate the photobleaching there
    for j=1:length(allresults.trial)
            traces = allresults.trial{j}.traces;
            trial_pb_ratios = photobleach_estimation(traces, 300); % Function used for original submission: trial_pb_ratios = photobleach_estimation(traces, 300);
            folder_pb_ratios = [folder_pb_ratios; trial_pb_ratios];
    end
    
    % Save the photobleaching ratios into the all results variable
    allresults.bleach = folder_pb_ratios;
    
    % Update/Save all results
    save([save_all_path fov_results(i).name], 'allresults');
end

% Remove Session Culture 29
delete([save_all_path '*29*.mat']);
