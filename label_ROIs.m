function [outputArg1,outputArg2] = label_ROIs(frame, ROIs)
    % Display average image
    figure;
    imagesc(frame);
    colormap(gray);

    % Show the ROIs on the average somArchon image
    for i=1:length(ROIs)
        B = bwboundaries(ROIs{i});
        hold on;
        visboundaries(B);
        hold on;
    
        % Use the mid rectangle of area function to find a center point
        [x, y] = centerRectPoly(ROIs{i});
    
        text(x-3, y, num2str(i), 'Color', 'yellow', 'FontSize', 14);
    end

end