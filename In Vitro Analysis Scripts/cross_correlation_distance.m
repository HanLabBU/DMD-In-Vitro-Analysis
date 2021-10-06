% This script will compute trace cross-correlation between neuron paris in
% a given field of view as a function of distance
% Author: Pierre Fabris

% Calculate cross-correlation for each neuron pair and then average
% their correlations across trials
%
% result
function [result] = cross_correlation_distance(trial_traces, ROIs)
    num_trials = size(trial_traces, 3);

    %% Saves all the cross correlation pair values here in each trial the useful
    % correlations are going to be in the upper triangle from the diagnonal.
    % Otherwise there would be redundant pair calculated
    cross_corr = zeros(length(ROIs), length(ROIs), num_trials);

    for trial=1:num_trials

        % Test that detrending an already detrended trace is ok
        traces = detrend(trial_traces(:, :, trial));
        cross_corr(:, :, trial) = corrcoef(traces);
    end

    % Average correlations across trials
    % Unique pairs are stored in the upper triangle of this matrix
    average_cross_corr = mean(cross_corr, 3);

    %% Calculate the centroids of each of the masks
    centroids = zeros(2, length(ROIs));
    for i=1:length(ROIs)
        boundary = bwboundaries(ROIs{i});
        x = boundary{1}(:, 2);
        y = boundary{1}(:, 1);

        % Remove duplicate vertices
        coors = [x y];
        coors = unique(coors, 'rows');
        
        polyin = polyshape(coors(:,1), coors(:,2));
        [centerx, centery] = centroid(polyin);

        % Save centroids
        centroids(:, i) = [centerx; centery];
    end
    
    % Save the centroids
    result.centroids = centroids;
    
    %% Calculate the distance between each neuron centroid pair
    % match the neuron pair with their cross-correlation
    % This variable is actually not used yet
    neuron_dist = zeros(length(centroids), length(centroids)); % Unique pairs are going to be the upper triangle of matrix

    % Store the distance and corresponding correlation
    dist_ave_corr = [];

    for i=1:length(centroids) - 1
        for j=i+1:length(centroids)
            x1 = centroids(1, i);
            y1 = centroids(2, i);
            x2 = centroids(1, j);
            y2 = centroids(2, j);

            dist = sqrt((x1 - x2).^2 + (y1 - y2).^2);
            neuron_dist(i, j) = dist;

            % Record the distance and the corresponding average correlation
            dist_ave_corr = [dist_ave_corr, [dist; average_cross_corr(i, j)]];
        
        end
    end

    % Save the cross correlation for each trial and the distance 
    for trial=1:num_trials
        result.trial{trial}.cross_correlation = cross_corr(:, :, trial);
    end
    result.neuron_dist = neuron_dist;
    
    % Save the average cross correlations
    result.average_cross_corr = average_cross_corr;
    
    % Perform linear regression
    dists = dist_ave_corr(1, :)';
    corrs = dist_ave_corr(2, :)';
    X = [ones(length(dists), 1) dists];
    b = X\corrs;
    regy = X*b;

    % Save the regressed values and the regressors
    result.best_fit_y = regy;
    result.regressors = b;
    
    % Calculate the R^2
    R2 = 1 - sum((corrs - regy).^2)/sum((corrs - mean(corrs)).^2);
    
    % Save the R^2 value
    result.R2 = R2;
    
    hold on;
    plot(dists, regy, '-', 'LineWidth', 3);
    message = ['slope=' num2str(b(2)) newline 'R^2 = ' num2str(R2)];
    annotation('textbox', [.2 .6 .9 .3], 'String', message, 'FitBoxToText', 'on');

end
