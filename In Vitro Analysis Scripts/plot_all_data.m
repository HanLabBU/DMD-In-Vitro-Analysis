% Used to plot all of the data

%% Setup for plotting the data from the allresults folder
% Folder to save all analysis data
save_all_path = '\\ad\eng\research\eng_research_handata\Pierre Fabris\DMD Project\All In Vitro Analysis\';


%% Generate SNR plots of each culture data
listings = dir(save_all_path);
results_files = {listings.name};
indi_files = results_files(contains(results_files, 'Individual'));

% Store all of the respective mask's SNR values as column vectors
indi_SNR = [];
wide_SNR = [];

% Match each individual mask file to its wide field
for i = 1:length(indi_files)
    indi_file = indi_files{i};
    fov_name = erase(indi_file, '_Individual Mask_allresults.mat');
    wide_file = results_files(contains(results_files, fov_name) & contains(results_files, 'Wide Field'));
    wide_file = wide_file{1};
    
    load([save_all_path indi_file]);
    indi_results = allresults;
    
    load([save_all_path wide_file]);
    wide_results = allresults;
    
    % Check for matching ROIs and store them 
    match_rois_idx = []; % Top row is individual mask; Bottom row is wide field
    for i=1:length(indi_results.roi_centroids)
        for j=1:length(wide_results.roi_centroids)
            x1 = indi_results.roi_centroids(1, i);
            y1 = indi_results.roi_centroids(2, i);
            x2 = wide_results.roi_centroids(1, j);
            y2 = wide_results.roi_centroids(2, j);
            
            pixdiff = sqrt((x1 - x2).^2 + (y1 - y2).^2);
            if pixdiff < 5
                match_rois_idx(:, :, end + 1) = [i; j];
                j;
                break;
            end
        end
    end
    match_rois_idx = match_rois_idx(:, :, 2:end);
        
    % Store the average SNR of all individual traces
     % TODO congregate all cutlure SNRs for individual
    for i=1:length(indi_results.snr.trial)
        trial_SNR = indi_results.snr.trial{i}.spike_snr(match_rois_idx(1, :));
        
        % Compute each neuron's average 
        ave_trial_SNR = [];
        for i=1:length(trial_SNR)
            ave_trial_SNR = [ave_trial_SNR; nanmean(trial_SNR{i})];
        end
        indi_SNR = horzcat_pad(indi_SNR,  ave_trial_SNR);
    end
    
    % Store the average SNR of all individual traces
    
    for i=1:length(wide_results.snr.trial)
        trial_SNR = wide_results.snr.trial{i}.spike_snr(match_rois_idx(2, :));
        
        % Compute each neuron's average 
        ave_trial_SNR = [];
        for i=1:length(trial_SNR)
            ave_trial_SNR = [ave_trial_SNR; nanmean(trial_SNR{i})];
        end
        
        wide_SNR = horzcat_pad(wide_SNR, ave_trial_SNR);
    end
    
    indiSNR = indi_SNR(:, end);
    wideSNR = wide_SNR(:, end);
    
    % Plot the bar plots with error bars
    figure('COlor','w'),,plot(indiSNR,'r'); hold on,plot(wideSNR,'k');
    legend indi wide;
    title([fov_name ' all trial SNR']);
    [h,p,ci,stats] = ttest(indiSNR,wideSNR);
    
    figure('COlor','w','Position', [ 300 300 300 350]);
    V1=nanmean(indiSNR);V1s=nanstd(indiSNR)./sqrt(length(indiSNR))
    V2=nanmean(wideSNR);V2s=nanstd(wideSNR)./sqrt(length(wideSNR))
    bar( [ 1 ], [V1],0.7,'FaceColor', [ 0.7 0.2 0.1]) , hold on,bar( [ 2 ], [V2],0.7,'FaceColor', [0.1 0.4 0.7]);
    set(gca,'Xtick', [ 1 2],'Xticklabel', {'DMD' ; 'Widefield'});
    errorbar([ 1 2], [ V1 V2], [V1s V2s],'.k','Linewidth', 2);
    axis tight;ylabel('Spike SNR');
    xlim([ 0.5 2.5]); ylim([3 5]);
    ax = gca;
    ax.XAxis.FontSize = 15;
    ax.YAxis.FontSize = 15;
    title([fov_name ' SNR Comparison(p= ' num2str(p) ')']);
    % Save the matched indices for the other illumination pattern
end

% Overall SNR between individual and wide field