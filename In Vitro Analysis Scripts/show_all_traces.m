
% % in vivo data
%cd('\\engnas.bu.edu\research\eng_research_handata\EricLowet\DMD\invivoDMD\')
%cd('\\engnas.bu.edu\research\eng_research_handata\EricLowet\DMD\comp_wide_indi\')
clc;
close all;
clear all;
 
% Use the scripts in the DMD scripts folder
addpath('.');

% in vitro data
cd('~/handata_server/Pierre Fabris/DMD Project/All In Vitro Analysis/');

% Folder to save figures
%save_fig_path = '\\engnas.bu.edu\research\eng_research_handata\Pierre Fabris\DMD Project\Data Fi226 - 67 (NaNs) = gures\';
save_fig_path = '~/handata_server/Pierre Fabris/DMD Project/Data Figures/';

% Get all of the analyzed data
ses=dir('*.mat');

% Find corresponding individual field
indiloc = [];

%Loop through every file and plot each trace
for file=1:length(ses)
    load(ses(file).name);

    % Loop through each trial
    for tr=1:length(allresults.trial)
        
        % Loop through each neuron
        for neuron=1:size(allresults.trial{tr}.traces, 2)
            trace = allresults.trial{tr}.traces(:, neuron);

            figure('Position', [0 0 2000 1500]);
            plot(trace);
            ylim([min(trace) - 40, max(trace) + 40]);
            title([ allresults.fov_name ' ' allresults.type ' Trial ' num2str(tr) ' Neuron ' num2str(neuron)]);
        end
    end
end
