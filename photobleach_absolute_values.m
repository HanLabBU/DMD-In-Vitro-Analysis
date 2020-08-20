function [result] = photobleach_absolute_values(traces, n)
	show_fig = 0;
    
    pb_diff = [];
    for i=1:size(traces, 2)
        trace = traces(:, i);
        
        % Determine the bounds for plotting the before and after filtered
        % traces
        min_y = min(trace) - 20;
        max_y = max(trace) + 20;
        
        % Show original trace
        if show_fig == 1
            figure();
            subplot(1, 2, 1);
            plot(trace);
            ylim([min_y max_y]);
        end
        
        trace = medfilt1(trace, 51);
            
        % Show median filtered trace
        if show_fig == 1
            subplot(1, 2, 2);
            plot(trace);
            ylim([min_y max_y]);
            sgtitle(['Neuron ' num2str(i)]);
        end
        init_mean_intensity = nanmean(trace(1:n)); %-767.7;
        last_mean_intensity = nanmean(trace(end-n:end)); %-767.7;
        
        
        pb_diff = [pb_diff, last_mean_intensity - init_mean_intensity];
    end
    
    result = pb_diff;
end