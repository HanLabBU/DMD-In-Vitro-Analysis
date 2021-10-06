% In vivo path
%cd('\\engnas.bu.edu\research\eng_research_handata\EricLowet\DMD\invivoDMD\')
%cd('\\engnas.bu.edu\research\eng_research_handata\EricLowet\DMD\comp_wide_indi\')

close all;
clear all;
clc;

addpath('.');

% In vitro path
cd('~/handata_server/Pierre Fabris/DMD Project/All In Vitro Analysis/');
%cd('D:\DMD Analysis Temp Data\');


% Scripts for needed functions
addpath('~/handata_server/EricLowet/DMD/main_analysis/');

ses=dir('*.mat');

% Destination folder for spike roasters
dest_dir = '~/handata_server/0Collaboration with Gabe Ocker/Invitro Voltage Data from Pierre/';

%%% search for wide
clear findwide
for id=1:length(ses) 
    if strfind(ses(id).name, 'Wide')>0
        findwide(id)=1;
    else
        findwide(id)=0;
    end
                                  
end
    
wideloc=find(findwide);

% Find corresponding individual field
indiloc = [];
for file=1:length(ses)
    if contains(ses(file).name, 'Individual') == 1
        load([ses(file).name]);
        indi_name = allresults.fov_name;

        for id=1:length(ses)
            if contains(ses(id).name, 'Wide') == 1 & contains(ses(id).name, indi_name) == 1           
                % Sanity check
                load([ses(id).name]);

                if strcmp(allresults.fov_name, indi_name) == 1
                    indiloc = [indiloc, file];
                end
            end
        end
    end
end


trial_size = 10000;
max_neurons = 0;
min_neurons = 14124124;
% Loop through each session
for id=1:length(wideloc)
    widefile=load(ses(wideloc(id)).name);
    indifile=load(ses(indiloc(id)).name);
    
    %Calculate min max neurons for each field of view
    if length(indifile.allresults.roi) > max_neurons
        max_neurons = length(indifile.allresults.roi);
    end
    
    if length(indifile.allresults.roi) < min_neurons
        min_neurons = length(indifile.allresults.roi);
    end
    
%     %Grab the DMD spike points
%     for i=1:length(indifile.allresults.trial)
%         spike_raster = indifile.allresults.roaster(:, (i*10000 - 10000 + 1):i*10000)';    
%         detrend_trace = indifile.allresults.orig_trace(:, (i*10000 - 10000 + 1):i*10000)';
%         trace_noise = indifile.allresults.trace_noise(:, i)';
%         
%         % Save the trials
%         trial_file = [dest_dir indifile.allresults.fov_name '/Individual Mask/raster_' num2str(i) '.mat'];
%         
%         % TODO temporary comment out
%         %save(trial_file, 'spike_raster', 'detrend_trace', 'trace_noise');
%     end
    
end
