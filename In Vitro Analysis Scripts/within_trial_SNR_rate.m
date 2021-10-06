
% In vivo path
%cd('\\engnas.bu.edu\research\eng_research_handata\EricLowet\DMD\invivoDMD\')
%cd('\\engnas.bu.edu\research\eng_research_handata\EricLowet\DMD\comp_wide_indi\')

close all;
clear all;

addpath('.');

% In vitro path
cd('~/handata_server/Pierre Fabris/DMD Project/All In Vitro Analysis/');
%cd('D:\DMD Analysis Temp Data\');

% Final figure folder
save_fig_path = '~/handata_server/Pierre Fabris/DMD Project/Data Figures/';

% Scripts for needed functions
addpath('~/handata_server/EricLowet/DMD/main_analysis/');

% Find the individual and wide field files
ses=dir('*.mat');
for id=1:length(ses) 
	if strfind(ses(id).name, 'Wide')>0
     		findwide(id)=1;
 	else
      		findwide(id)=0;
 	end
end

wideloc=find(findwide);
indiloc=find(~findwide);

if length(wideloc) ~= length(indiloc)
	disp('Number of individual and wide field FOVs are not equal');
	pause;
end

% Read each FOV pair files
for id=1:length(indiloc)
	indifile = load(ses(indiloc(id)).name)
	widefile = load(ses(wideloc(id)).name)
	


end


