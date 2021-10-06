% This function takes in a mask with a single area and finds the extreme
% values in each dimension. Then calculate the midpoint for each dimension
function [x, y] = centerRectPoly(mask)
    rows = size(mask, 1);
    cols = size(mask, 2);
    
    min_r = rows; max_r = 0;
    min_c = cols; max_c = 0;
    
    % Loop through entire mask and find which entries have a logical true
    % value
    for i = 1:rows
        for j=1:cols
            if mask(i, j) == 1
                
                %Update the extremes
                if i < min_r min_r = i; end
                if i > max_r max_r = i; end
                if j < min_c min_c = j; end
                if j > max_c max_c = j; end
                
            end
        end
    end
    
    % Find the midpoint of this rectangle
    x = (min_c + max_c)/2;
    y = (min_r + max_r)/2;
end