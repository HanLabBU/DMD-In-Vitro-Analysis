function [outputArg1,outputArg2] = label_ROIs(frame, ROIs)
    % Display average image
    figure('Position', [0 0 2000 1500]);
    imagesc(frame);
    colormap(gray);

    % Show the ROIs on the average somArchon image
    for i=1:length(ROIs)
        B = bwboundaries(ROIs{i});
        hold on;
        visboundaries(B, 'LineStyle', ':', 'EnhanceVisibility', false, 'Color', 'yellow');
        hold on;
    
        % Use the mid rectangle of area function to find a center point
        [x, y] = centerRectPoly(ROIs{i});
    
        text(x-5, y, num2str(i), 'Color', 'cyan', 'FontSize', 18);
    end

end