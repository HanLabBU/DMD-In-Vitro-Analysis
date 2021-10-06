% Takes trial traces and calculates the percent change in
% intensity between the first n points and the n 100 points

% Author: Pierre Fabris

% Arguments: 
% traces - all of the traces in a given trial
% n - the number of points at the beginning or end of trace to calculate
% the intensity ratio
%
% Returns an array of all the photobleach ratios
function [result] = photobleach_estimation_final(traces, n)
    photo_bleach_ratio = [];
  noise_floor= 767.7;
  time_span=5000;
    % Loop through each neuron
    for i=1:size(traces, 2)
            trace = traces(:, i);
            
            % Apply median filter to filter out spikes
            % This filter is not changing much of the overall plots
            trace = medfilt1(trace, 51);
            if ~any(isnan(trace))
                if length(traces)<= time_span
            init_mean_intensity = nanmean(trace(1:n)) - noise_floor;
            last_mean_intensity = nanmean(trace(end-n:end)) - noise_floor;
                else
                init_mean_intensity = nanmean(trace(1:n)) - noise_floor;
            last_mean_intensity = nanmean(trace( time_span-n: time_span)) - noise_floor;     
                end
            photo_bleach_ratio = [photo_bleach_ratio; [last_mean_intensity/init_mean_intensity]];
            
            
            end
            
    end
    
    result = photo_bleach_ratio;
    
end
