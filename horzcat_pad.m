function [final_array] = horzcat_pad(main_array, new_array)
    % New folder has a boatload of traces
    if size(new_array, 1) > size(main_array, 1)
        pad = NaN(length(new_array) - size(main_array, 1), size(main_array, 2));
        main_array = [ [main_array ; pad], new_array];
    
    % The new folder has fewer traces
    elseif size(main_array, 1) > length(new_array)
        pad = NaN(size(main_array, 1) - length(new_array), size(new_array, 2));
        new_array = [new_array ; pad];
        main_array = [main_array, new_array];
        
        % The number of traces was miraculously the same
    else
        main_array = [main_array, new_array];
    end
    
    final_array = main_array;
end

