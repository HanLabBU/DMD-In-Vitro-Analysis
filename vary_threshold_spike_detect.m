clear all;
close all;
clc;

%% Author notes for using script
% Will be using the trial{1-3}.traces to do the spike detection


% Save figure directory
save_fig_path = '~/handata_server/Pierre Fabris/DMD Project/Data Figures/';


% Data directory
invitro_data_path = '/home/pierfier/Projects/DMD Local Data/';


ses=dir([invitro_data_path '/*.mat']);


% Find all of the 
% Search for wide
clear findwide
for id=1:length(ses) 
 if   strfind(ses(id).name, 'Wide')>0
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
        indifile = load([ses(file).name]);
        indi_name = indifile.allresults.fov_name;

        for id=1:length(ses)
            if contains(ses(id).name, 'Wide') == 1 & contains(ses(id).name, indi_name) == 1           
                % Sanity check
                widefile = load([ses(id).name]);
                if strcmp(widefile.allresults.fov_name, indi_name) == 1 & strcmp(widefile.allresults.type, 'wide') == 1
                    indiloc = [indiloc, file];
                end
            end
        end
    end
end

% Sanity check that matching pairs are found
if length(indiloc) ~= length(wideloc)
	disp('Not all of the data found are matched pairs');
	pause;
end

% Store the SNRs with the corresponding threshold, so {<threshold>, <SNRs>}
indi_SNR = {};
wide_SNR = {};

% Range of threshold values to test
up_thres = [0.1:.1:8];

% Loop through all of the data files
for id=1:length(indiloc)
	indifile = load(ses(indiloc(id)).name);
	widefile = load(ses(wideloc(id)).name);
	indiresults = indifile.allresults;
	wideresults = widefile.allresults;

	indi_temp = [];
	% Calculate individual DMD SNRS 
	parfor trial = 1:length(indiresults.trial)
		
	end	

	% Calculate wide field SNRS 
	parfor trial = 1:length(wideresults.trial)
		
	end	

end
