function [outputArg1,outputArg2] = display_cells_in_frame(frame, ROIs, scale_bar)
    % Display average image
    figure('Renderer', 'painters', 'Position', [0 0 2000 1500]);
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
    
        text(x-40, y, num2str(i), 'Color', 'cyan', 'FontSize', 18);
    end
    
    % Plot the scale bar
    plot([scale_bar(1) scale_bar(1)+scale_bar(3)], [scale_bar(2) scale_bar(2)], '-w', 'LineWidth', scale_bar(4));
end
