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
    
        %text(x-40, y, num2str(i), 'Color', 'cyan', 'FontSize', 18);
    end
    
    % Plot the scale bar
    sc=0.1625*2; % micrometer per pixel 40x
    ten_um_pixels = 10./sc;
    posx = 360;
    plot([posx posx+ten_um_pixels], [500 500], '-w', 'LineWidth', 4);
    
end
