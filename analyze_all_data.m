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

%% TODO make a structure that pairs all of the data?
% Also, save the ROIs, the original traces, the name of the FOV folder, and
% the mask type
for i = 1:length(FOV_dirs)
    fov = FOV_dirs{i};
    listings = dir(fov);
    files = {listings.name};
   
    
    
    % Check if all results file exists for this FOV in the saved folder
    allresults_file = [strrep(erase(fov, root_path), '\', '_') 'allresults'];
    allresults.fov_name = strrep(erase(fov, {root_path, 'Individual Mask', 'Wide Field'}), '\', '');
    
    % Check if all results already exists and load its contents
    % After this point all of the allresults files from root_path should
    % be created
    if exist([save_all_path allresults_file '.mat'])
        load([save_all_path allresults_file]);
    end
    
    % Read in all trials and save original traces
    trace_files = files(contains(files, 'traces'));
    allresults.trial = {};
    
    for i=1:length(trace_files)
        load([fov trace_files{i}]);
        allresults.trial{i}.traces = traces;
    end
    
    % Save all the ROIs
    roi_file = files(contains(files, 'ROIs'));
    load([fov roi_file{1}]);
    allresults.rois = ROIs;
    
    % Save the type of mask
    if contains(fov, 'Individual Mask')
        allresults.type = 'indi';
    elseif contains(fov, 'Wide Field')
        allresults.type = 'wide';
    end
    
    save([save_all_path allresults_file], 'allresults');
end

%% Calculate the cross correlations and distance between neurons
for i = 1:length(FOV_dirs)
    fov = FOV_dirs{i};
    listings = dir(fov);
    files = {listings.name};
    
    % Check if allresults file exists for this FOV in the saved folder
    allresults_file = [strrep(erase(fov, root_path), '\', '_') 'allresults']
    load([save_all_path allresults_file]);
    
    
    % Compile all of the traces in the trial
    max_trial = 0;
    for i=1:length(allresults.trial)
        max_trial = max(max_trial, size(allresults.trial{i}.traces, 1));
    end
    
    trial_traces = NaN(max_trial, size(allresults.trial{1}.traces, 2));
    for i=1:length(allresults.trial)
        trial_traces(1:size(allresults.trial{i}.traces, 1), :, i) = allresults.trial{i}.traces;
    end
    
    % TODO will need to make sure all of the results data is saved
    cross_corr_result = cross_correlation_distance(trial_traces, allresults.ROIs);
    
    allresults.roi_centroids = cross_corr_result.centroids;
    allresults.neuron_dist = cross_corr_result.neuron_dist;
    
    save([save_all_path allresults_file], 'allresults');
end

%% Detect spikes and calculate SNR
% At the moment, different tweaks are going to be made to the SNR algorithm
for i = 1:length(FOV_dirs)
    fov = FOV_dirs{i};
    listings = dir(fov);
    files = {listings.name};
   
    % Check if all results file exists for this FOV in the saved folder
    allresults_file = [strrep(erase(fov, root_path), '\', '_') 'allresults'];
    load([save_all_path allresults_file]);
    
    
    for i=1:length(allresults.trial)
        traces = allresults.trial{i}.traces;
        
        % Perform the spike_detection_SNR
        result = spike_detect_SNR_v3(traces);
        
        % Save all of the SNR results data into the allresults structure
        allresults.orig_trace{:, i} = result.orig_trace;
        allresults.denoise_trace{:, i} = result.denoise_trace;
        allresults.trace_ws{:, i} = result.trace_ws;
        allresults.orig_traceDN{:, i} = result.orig_traceDN;
        allresults.roaster{:, i} = result.roaster;
        allresults.spike_snr{:, i} = result.spike_snr;
        allresults.spike_amplitude{:, i} = result.spike_amplitude;
        allresults.spike_idx{:, i} = result.spike_idx;
        allresults.trace_noise{:, i} = result.trace_noise;        
    end
    
    % Save all the SNR results 
    save([save_all_path allresults_file], 'allresults');
end

%% Calculate the photobleaching ratios for each FOV folder
% Store all of the photobleaching ratios
all_pb_ratios = [];
pb_plot_labels = {};
for i = 1:length(FOV_dirs)
    fov = FOV_dirs{i};
    listings = dir(fov);
    files = {listings.name};
    isdirs = [listings.isdir];
    
    % Check if all results file exists for this FOV in the saved folder
    allresults_file = [strrep(erase(fov, root_path), '\', '_') 'allresults'];
    load([save_all_path allresults_file])
    
    % Store all of the ratios in a single folder
    folder_pb_ratios = [];
    
    % Trial mean/std of photobleaching ratio
    trial_pb_averages_std = [];
    
    % Load each trial's traces
    for i=1:length(files)
        % Make sure that it is the traces mat file
        if contains(files{i}, 'traces') 
            load([fov files{i}]);
            trial_pb_ratios = photobleach_estimation(traces, 300);
            folder_pb_ratios = [folder_pb_ratios; trial_pb_ratios];
            
            trial_pb_averages_std = [trial_pb_averages_std, [nanmean(trial_pb_ratios); std(trial_pb_ratios)]];
        end
    end
    
    % Average and std of pb ratio for entire folder
    average_folder_pb = [nanmean(folder_pb_ratios); std(folder_pb_ratios)];
    
    % Save the photobleaching ratios into the all results variable
    allresults.bleach = folder_pb_ratios;
    
    % Add folder to all ratios with averaging across trials
    all_pb_ratios = horzcat_pad(all_pb_ratios, nanmean(folder_pb_ratios, 2));
    pb_plot_labels{end + 1} = [allresults.fov_name ' ' allresults.type];
    
    % Update/Save all results
    save([save_all_path allresults_file], 'allresults');
    
end


% Plot all of the folder photobleachings
figure;
boxplot(all_pb_ratios, pb_plot_labels);
ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
xtickangle(ax, 45);
ylabel('(initial intensity)/(final intensity)');


