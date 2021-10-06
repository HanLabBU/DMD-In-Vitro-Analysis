function hippocampus_event_detection(traces, up_threshold, down_threshold)

    if nargin<3 || isempty(down_threshold) 
        down_threshold = 3;
    end
    
    if nargin<2 || isempty(up_threshold) 
        up_threshold = 4;
    end
    
    % Camera Frequency
    FS = 500; % Hz
    
    show_figure = 1;
    show_snr_figure = 0;
    refine_event_with_snr_cluster = 1;

    event_parameter.moving_window = 201;
    event_parameter.pre_peak_data_point = 1;
    event_parameter.post_peak_data_point = 1;
    event_parameter.event_moving_window = 201; % data points
       
    event_parameter.noise_threshold = 5; %7
    event_parameter.noise_pre_extension = 0; % data points
    event_parameter.noise_post_extension = 3; % data points
    event_parameter.noise_extension = 3; % data points
    event_parameter.noise_moving_window = 11; % data points
    
    event_parameter.snr_threshold = 0;
%     event_parameter.snr_window = 100;
    event_parameter.refine_threshold = 4; % standard deviation
    event_parameter.down_threshold = down_threshold;
    event_parameter.up_threshold = up_threshold;

    whole_tic = tic;
        
        result.denoise_traces = [];
        result.traces = traces;
        
        for roi_idx=1:size(traces,2)
            clear event;
            event.idx=[];
            event.amplitude=[];
            event.snr=[];
            current_trace = traces(:,roi_idx);
                       
            trace_time = 1:size(traces(:, 1));
            
            f_trace = eegfilt(current_trace', FS,5,0)';
            u_f_trace = get_upper_trace(f_trace,event_parameter.moving_window);
            
            d_u_f_trace = diff(u_f_trace);
            d_u_f_trace = [0;d_u_f_trace];

            noise_idx_list = find_noise_idx(f_trace,event_parameter.moving_window,event_parameter.noise_moving_window,event_parameter.noise_threshold,event_parameter.noise_extension);
            result.roi(roi_idx).noise_idx = noise_idx_list;
            
            denoise_trace = current_trace;
            for idx=1:numel(noise_idx_list)
                current_noise_idx = noise_idx_list(idx);
                if current_noise_idx-event_parameter.noise_pre_extension>0
                    denoise_trace(current_noise_idx-event_parameter.noise_pre_extension:min(current_noise_idx+event_parameter.noise_post_extension,length(current_trace))) = nan;
                end
            end
            result.denoise_traces = cat(2,result.denoise_traces,denoise_trace);
            
%             d_denoise_trace = diff(denoise_trace);
%             d_denoise_trace = [0;d_denoise_trace];

            l_f_trace = get_lower_trace(f_trace,event_parameter.moving_window);
            d_l_f_trace = diff(l_f_trace);
            d_l_f_trace = [0;d_l_f_trace];
            
            current_trace_noise = 2*std(l_f_trace); 
            
            event_parameter.up_threshold_value = event_parameter.up_threshold*nanstd(d_l_f_trace);
            event_parameter.down_threshold_value = event_parameter.down_threshold*nanstd(d_l_f_trace);
            
            event.event_parameter.up_threshold_value = event_parameter.up_threshold_value;           
            event.event_parameter.down_threshold_value = event_parameter.down_threshold_value;
            
            event.trace = current_trace;
            event.denoise_trace = denoise_trace;
            
            pre_d_trace = d_u_f_trace;
            if event_parameter.pre_peak_data_point>0
                for idx=1:event_parameter.pre_peak_data_point
                    shifted_d_trace = [zeros(idx,1);d_u_f_trace(1:end-idx)];
                    shifted_d_trace(shifted_d_trace<0) = 0;
                    pre_d_trace = pre_d_trace+shifted_d_trace;
                end
            end
            
            post_d_trace = d_u_f_trace;
            if event_parameter.post_peak_data_point>0
                for idx=1:event_parameter.post_peak_data_point
                    shifted_d_trace = [d_u_f_trace(idx+1:end);zeros(idx,1)];
                    shifted_d_trace(shifted_d_trace>0) = 0;
                    post_d_trace = post_d_trace+shifted_d_trace;
                end
            end
            
            up_idx_list = find(pre_d_trace>(nanmean(d_u_f_trace)+event_parameter.up_threshold_value));
            
            for up_idx=up_idx_list'
                if up_idx>2 & (up_idx+1)<=numel(d_u_f_trace) & d_u_f_trace(up_idx)>0 & d_u_f_trace(up_idx+1)<0 & post_d_trace(up_idx+1)<(nanmean(d_u_f_trace)-event_parameter.down_threshold_value) & ~isnan(denoise_trace(up_idx))
           
                    % SNR
                    peak_intensity = current_trace(up_idx);
                    pre_peak_intensity_1 = current_trace(up_idx-1);
                    pre_peak_intensity_2 = current_trace(up_idx-2);
                    
%                     snr_window_start = max(1,up_idx-event_parameter.snr_window);
%                     snr_window_end = min(numel(l_f_trace),up_idx+event_parameter.snr_window);
%                     current_trace_noise = 2*std(l_f_trace(snr_window_start:snr_window_end));

                    current_signal_intensity = max(peak_intensity-pre_peak_intensity_1,peak_intensity-pre_peak_intensity_2);
                    current_snr = current_signal_intensity/current_trace_noise;
                    
                    if current_snr>=event_parameter.snr_threshold
                        event.idx = cat(1,event.idx,up_idx);
                        event.amplitude = cat(1,event.amplitude,current_signal_intensity);
                        event.snr = cat(1,event.snr,current_snr);
                    end
                    
                end
            end

            event.time = trace_time(event.idx);
            event.roaster = zeros(size(current_trace));
            event.roaster(event.idx) = 1;
            event.trace_noise = current_trace_noise;
            event.snr_threshold = event_parameter.snr_threshold;
            
            [event.bursting_idx,event.bursting_event_idx,event.bursting_window_idx,event.bursting_window_time,event.bursting_duration] = find_bursting(event.idx,event.time);
            event.bursting_time = trace_time(event.bursting_idx);

            result.roi(roi_idx).event = event;
            
            
            if refine_event_with_snr_cluster==1
                current_all_snr = event.snr;
                if numel(current_all_snr)>1
                    idx = kmeans(current_all_snr,2,'start',[mean(current_all_snr)-2*std(current_all_snr);mean(current_all_snr)+2*std(current_all_snr)]);
                    if sum(idx==1)>sum(idx==2)
                        major_cluster_idx = find(idx==1);
                    elseif sum(idx==1)<sum(idx==2)
                        major_cluster_idx = find(idx==2);
                    else
                        major_cluster_idx = 1:numel(current_all_snr);
                    end
                
                    refine_snr_threshold = mean(current_all_snr(major_cluster_idx))-event_parameter.refine_threshold*std(current_all_snr(major_cluster_idx));
                else
                    refine_snr_threshold = 0;
                end
                
                refine_idx = find(current_all_snr>=refine_snr_threshold);
                refined_event.idx = event.idx(refine_idx);
                refined_event.amplitude = event.amplitude(refine_idx);
                refined_event.snr = event.snr(refine_idx);
                refined_event.time = trace_time(refined_event.idx);
                refined_event.roaster = zeros(size(current_trace));
                refined_event.roaster(refined_event.idx) = 1;
                refined_event.trace_noise = current_trace_noise;
                refined_event.snr_threshold = refine_snr_threshold;

                [refined_event.bursting_idx,refined_event.bursting_event_idx,refined_event.bursting_window_idx,refined_event.bursting_window_time,refined_event.bursting_duration] = find_bursting(refined_event.idx,refined_event.time);
                refined_event.bursting_time = trace_time(refined_event.bursting_idx);

                result.roi(roi_idx).refined_event = refined_event;
                
                remove_idx = find(current_all_snr<refine_snr_threshold);
                removed_event.idx = event.idx(remove_idx);
                removed_event.amplitude = event.amplitude(remove_idx);
                removed_event.snr = event.snr(remove_idx);
                removed_event.time = trace_time(removed_event.idx);
                removed_event.roaster = zeros(size(current_trace));
                removed_event.roaster(removed_event.idx) = 1;
                removed_event.trace_noise = current_trace_noise;
                removed_event.snr_threshold = refine_snr_threshold;

                result.roi(roi_idx).removed_event = removed_event;
                    
            end


        end
        
        result.event_parameter = event_parameter;

    if show_figure==1
        figure;
        x_axis = gca;
        for roi_idx=1:size(traces,2)
            current_trace = result.roi(roi_idx).event.trace;
            current_denoise_trace = result.roi(roi_idx).event.denoise_trace;
            
            
            current_m_trace = movmean(current_trace,event_parameter.moving_window);
            current_l_trace = get_lower_trace(current_trace,event_parameter.moving_window);
            current_std_trace = movstd(current_l_trace,event_parameter.moving_window);

            
%             subplot(2,1,1)
            hold on
            plot(x_axis,current_trace,'color',[0.7,0.7,0.7])
            plot(x_axis,current_denoise_trace,'color','b')
            plot(x_axis,current_m_trace,'color','k')
            rectangle_heigh = (max(current_trace)-min(current_trace));
            bursting_time = 0;
            if refine_event_with_snr_cluster==1
                current_refined_event_idx = result.roi(roi_idx).refined_event.idx;
                current_removed_event_idx = result.roi(roi_idx).removed_event.idx;
                scatter(x_axis(current_refined_event_idx),current_denoise_trace(current_refined_event_idx),8,'r','filled');
                scatter(x_axis(current_removed_event_idx),current_denoise_trace(current_removed_event_idx),8,'g','filled');
   
                if ~isempty(result.roi(roi_idx).refined_event.bursting_window_idx)
                    bursting_count = size(result.roi(roi_idx).refined_event.bursting_window_idx,1);
                    for idx = 1:bursting_count
                        current_bursting_start = x_axis(result.roi(roi_idx).refined_event.bursting_window_idx(idx,1));
                        current_bursting_end = x_axis(result.roi(roi_idx).refined_event.bursting_window_idx(idx,2));
                        rectangle_width = current_bursting_end-current_bursting_start;
                        bursting_time = bursting_time+rectangle_width;
                        rectangle('Position',[current_bursting_start min(current_trace) rectangle_width rectangle_heigh],'FaceColor',[rand(1) rand(1) rand(1) 0.3],'EdgeColor','none')
                    end
                else
                    bursting_count = 0;
                end
                
                avg_refined_snr = mean(result.roi(roi_idx).refined_event.snr);
                avg_removed_snr = mean(result.roi(roi_idx).removed_event.snr);
                title({['ROI ',num2str(roi_idx)],['SNR: ',num2str(avg_refined_snr),' / ',num2str(avg_removed_snr)],['Bursting: ',num2str(bursting_count),' / ',num2str(bursting_time)]})
            else
                current_event_idx = result.roi(roi_idx).event.idx;
                scatter(x_axis(current_event_idx),current_denoise_trace(current_event_idx),8,'r','filled');
                
                if ~isempty(result.roi(roi_idx).event.bursting_window_idx)
                    bursting_count = size(result.roi(roi_idx).event.bursting_window_idx,1);
                    for idx = 1:size(result.roi(roi_idx).event.bursting_window_idx)
                        current_bursting_start = x_axis(result.roi(roi_idx).event.bursting_window_idx(idx,1));
                        current_bursting_end = x_axis(result.roi(roi_idx).event.bursting_window_idx(idx,2));
                        rectangle_width = current_bursting_end-current_bursting_start;
                        bursting_time = bursting_time+rectangle_width;
                        rectangle('Position',[current_bursting_start min(current_trace) rectangle_width rectangle_heigh],'FaceColor',[rand(1) rand(1) rand(1) 0.3],'EdgeColor','none')
                    end
                else
                    bursting_count = 0;
                end
                
                avg_snr = mean(result.roi(roi_idx).event.snr);
                title({['ROI ',num2str(roi_idx)],['SNR: ',num2str(avg_snr)],['Bursting: ',num2str(bursting_count),' / ',num2str(bursting_time)]})
            end
            xlim([x_axis(1) x_axis(end)])
            xlabel('Time (sec)')
            
%             subplot(2,1,2)
%             hold on
%             plot(x_axis,all_isi,'b')
%             scatter(event_time(bursting_event_idx),zeros(numel(bursting_event_idx),1),8,'r','filled')
%             xlim([x_axis(1) x_axis(end)])
%             xlabel('Time (sec)')
%             ylabel('ISI (sec)')
            
            if show_snr_figure==1
                all_snr = zeros(size(current_trace));
                all_snr(result.roi(roi_idx).event.idx) = result.roi(roi_idx).event.snr;
                figure
                subplot(2,1,1)
                hold on
                plot(x_axis,all_snr,'color','b')
                if refine_event_with_snr_cluster==1
                    yline(result.roi(roi_idx).refined_event.snr_threshold,'r');
                end
                xlim([x_axis(1) x_axis(end)])
                xlabel('Time (sec)')
                ylabel('SNR')
                hold off                
                subplot(2,1,2)
                hist(result.roi(roi_idx).event.snr)
                if refine_event_with_snr_cluster==1
                    xline(result.roi(roi_idx).refined_event.snr_threshold,'r');
                end
                xlabel('SNR')
                ylabel('Event count')
                
                sgtitle(['ROI ',num2str(roi_idx)])
            end
            
            
        end
    end
    
    if ~isempty(error_file)
        fprintf('----------\n')
        fprintf('Error:\n')
        for idx=1:size(error_file,1)
            fprintf('\t%s\n',error_file{idx});          
        end
        fprintf('----------\n')
    end
    
    fprintf(['Total time: ',num2str(toc(whole_tic)),' seconds.\n']);

end

function noise_idx_list = find_noise_idx(current_trace,trace_moving_window,noise_moving_window,noise_threshold,noise_extension)

    m_trace = movmean(current_trace,trace_moving_window);
    lower_current_trace = current_trace;
    % replace the part above moving average with moving average
    idx = find(lower_current_trace>m_trace);
    lower_current_trace(idx)=m_trace(idx);

    movstd_lower_current_trace = movstd(lower_current_trace,noise_moving_window);
    noise_idx_list = find(isoutlier(movstd_lower_current_trace,'gesd')==1);
    %noise_idx_list = find(movstd_lower_current_trace>(mean(movstd_lower_current_trace)+noise_threshold*std(movstd_lower_current_trace)));

    % connect noise index
    noise_idx_list = sort(noise_idx_list);
    d_noise_idx_list = diff(noise_idx_list);
    noise_extension_idx = find(d_noise_idx_list>1 & d_noise_idx_list<noise_extension);
    if ~isempty(noise_extension_idx)
        for idx=1:numel(noise_extension_idx)
            current_idx = noise_extension_idx(idx);
            noise_idx_list = cat(1,noise_idx_list,[noise_idx_list(current_idx):noise_idx_list(current_idx+1)]');
        end
    end

    noise_idx_list = unique(noise_idx_list);

end

function [bursting_idx,bursting_event_idx,bursting_window_idx,bursting_window_time,bursting_duration] = find_bursting(event_idx,event_time)

    bursting_threshold = 0.03; % second
    event_isi = diff(event_time);
    bursting_event_idx = find(event_isi<bursting_threshold)+1;
    bursting_event_idx = unique([bursting_event_idx;bursting_event_idx-1]);
    bursting_idx = event_idx(bursting_event_idx);
    
    bursting_roaster = zeros(size(event_idx));
    bursting_roaster(find(event_isi<bursting_threshold)+1) = 1;
    bursting_start = find(diff(bursting_roaster)==1);
    bursting_end = find(diff(bursting_roaster)==-1);
    
    if ~isempty(bursting_start) & ~isempty(bursting_end)
        if bursting_start(1)>bursting_end(1)
            bursting_start = [1;bursting_start];
        end

        if bursting_start(end)>bursting_end(end)
            bursting_end = [bursting_end;numel(event_idx)];
        end

        bursting_window_idx = [event_idx(bursting_start) event_idx(bursting_end)];
        bursting_window_time = [event_time(bursting_start) event_time(bursting_end)];
        bursting_duration = bursting_window_time(:,2)-bursting_window_time(:,1);
    elseif isempty(bursting_start) & ~isempty(bursting_end)
        bursting_window_idx = [event_idx(1) event_idx(bursting_end)];
        bursting_window_time = [event_time(1) event_time(bursting_end)];
        bursting_duration = bursting_window_time(:,2)-bursting_window_time(:,1);
    elseif ~isempty(bursting_start) & isempty(bursting_end)
        bursting_window_idx = [event_idx(bursting_start) event_idx(numel(event_idx))];
        bursting_window_time = [event_time(bursting_start) event_time(numel(event_idx))];
        bursting_duration = bursting_window_time(:,2)-bursting_window_time(:,1);
    else
        bursting_window_idx = [];
        bursting_window_time = [];
        bursting_duration = [];
    end

end

function upper_trace = get_upper_trace(current_trace,trace_moving_window)

    m_trace = movmean(current_trace,trace_moving_window);
    upper_trace = current_trace;
    % replace the part below moving average with moving average
    idx = find(upper_trace<m_trace);
    upper_trace(idx)=m_trace(idx);
end

function lower_trace = get_lower_trace(current_trace,trace_moving_window)

    m_trace = movmean(current_trace,trace_moving_window);
    lower_trace = current_trace;
    % replace the part below moving average with moving average
    idx = find(lower_trace>m_trace);
    lower_trace(idx)=m_trace(idx);
end