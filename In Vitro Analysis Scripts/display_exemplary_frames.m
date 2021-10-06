% -- This script is used to get the image frames and then overlay cell boundaries

%% Setup shared variables
addpath("~/handata_server/EricLowet/DMD/preprocess/");


culture_7_videos_path = '~/handata_server/eng_research_handata3/DMD data/Good/Culture individual vs large fov/Culture 7/';

culture_analysis_data_path = '~/handata_server/Pierre Fabris/DMD Project/All In Vitro Analysis/';

figure_path = '~/handata_server/Pierre Fabris/DMD Project/Data Figures/Exemplary Trial or Sessions/SVG Format/';


sc= 362/1152; % Old value was 40x 0.1625*2; % micrometer per pixel 40x

%% Read and plot the targeted Culture 7 FOV


% This read function can only work for the DMD camera
% Also use only the first 100 frames
targeted_frame = read_dcimg_eric([culture_7_videos_path 'FoV_00001.dcimg'], 100);
targeted_frame = mean(targeted_frame, 3);

% Save the targeted ROIs
culture_7_indi_data = load([culture_analysis_data_path 'Culture 7_Individual Mask_allresults.mat']);
targeted_ROIs = culture_7_indi_data.allresults.roi;

% 20um scale bar 
scale_bar = [360, 500, 20./sc, 4]; % Array is x, y, length, then width
display_cells_in_frame(targeted_frame, targeted_ROIs, scale_bar);

saveas(gcf, [figure_path 'targeted_frame_cell_boundary.svg']);


%% Read and plot the widefield Culture 7 FOV
widefield_frame = read_dcimg_eric([culture_7_videos_path 'FoV_00002.dcimg'], 100);
widefield_frame = mean(widefield_frame, 3);

% Save the widefield ROIs
culture_7_wide_data = load([culture_analysis_data_path 'Culture 7_Wide Field_allresults.mat']);
widefield_ROIs = culture_7_wide_data.allresults.roi;

display_cells_in_frame(widefield_frame, widefield_ROIs, scale_bar);

saveas(gcf, [figure_path 'widefield_frame_cell_boundary.svg']);


%% Read and plot the GFP image of Culture 7 FOV
data = load([culture_7_videos_path 'FoV_001.mat']);;
gfp_frame = data.originalImage;

gfp_ROIs = {0, 0};

% 20um scale bar 
scale_bar = [360, 500, 20./sc, 4]; % Array is x, y, length, then width
display_cells_in_frame(gfp_frame, gfp_ROIs, scale_bar);

saveas(gcf, [figure_path 'gfp_frame_cell_boundary.svg']);


