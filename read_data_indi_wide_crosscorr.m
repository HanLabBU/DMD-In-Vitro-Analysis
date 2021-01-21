
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

ses=dir('*.mat');

%%% search for wide
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

% Original code for finding the corresponding wide field mask
%wideloc=find(findwide);
%%%% search for fitting indi
%indiloc=[];
%for id2=1:length(wideloc)
%widename=ses(wideloc(id2)).name;
%    for id=1:length(ses) 
%        if id~=wideloc(id2)
%            indiname=   ses(id).name;
%            vx2= strfind(widename,'wide');;vx= strfind(indiname,'fov');
%            if length(indiname) > vx2(1)-1
%                if length(strfind(indiname(1:vx2(1)-1), widename(1:vx2(1)-1)))>0
%                    indiloc=[indiloc, id];
%                end
%            end
%        end
%    end
%end

 sc= 0.1625*2 %362/1152;% Value if it was 40x for whole system 0.1625*2; % micrometer per pixel 40x
%% seelct matching ROI
%indiB=[];wideB=[]; indiSNR=[]; wideSNR=[];
cindi=[];cwide=[]; rdist=[];
cindiS=[];cwideS=[];  allCx=[];allCm=[];
for id=1:length(wideloc)
    widefile=load(ses(wideloc(id)).name);
    indifile=load(ses(indiloc(id)).name);
    clear wrm % ROI centroid
    for id2= 1:length(widefile.allresults.roi)
        [ x y]=find(widefile.allresults.roi{id2});
        wrm(:,id2)= round(mean([x , y]));end
    clear irm
    for id2= 1:length(indifile.allresults.roi)
        [ x y]=find(indifile.allresults.roi{id2});
        irm(:,id2)= round(mean([x , y]));end  
    mROI=[];
    for id3=1:size(irm,2) % matching ROI
        cents=irm(:,id3);
        pxdiff=(sqrt(sum(bsxfun(@minus, wrm, cents).^2)));
    
      wloc=find(pxdiff<8);
      if ~isempty(wloc)
    mROI=[mROI; [ id3 wloc]];end
    end
    
    %%% Put in matrix%%%
 %subthresIndi= indifile.allresults.trace_ws;
  % subthresWide= widefile.allresults.trace_ws;  
    subthresIndi= indifile.allresults.orig_traceDN;
   subthresWide= widefile.allresults.orig_traceDN;

   % subthresIndi= subthresIndi-fastsmooth( subthresIndi,10,1,1);
    %  subthresWide= subthresWide-fastsmooth( subthresWide,10,1,1);
nsel=mROI(:,1);
clear allCmax  ROIdist
for ind1=1:length(nsel)
for ind2=1:length(nsel)
%  A1=  zscore(subthresIndi(nsel(ind1),:));A2=  zscore(subthresIndi(nsel(ind2),:));
   % [c,lags]= xcorr(A1,A2,10,'Coeff');
        A1=  (subthresIndi(nsel(ind1),:));A2=  (subthresIndi(nsel(ind2),:));
    [CC]=corrcoef(A1,A2, 'Rows', 'complete');
    allCmax(ind1,ind2)= CC(1,2);%max(abs(c));
   ROIdist(ind1,ind2)= sqrt(sum((irm(:,nsel(ind1)).*sc-irm(:,nsel(ind2)).*sc).^2));
end;end
    allCmax(allCmax>0.99)=NaN;
    
clear allCmaxW  
for ind1=1:length(nsel)
for ind2=1:length(nsel)
%  A1=  zscore(subthresWide(nsel(ind1),:));A2=  zscore(subthresWide(nsel(ind2),:));
    A1=  (subthresWide(nsel(ind1),:));A2=  (subthresWide(nsel(ind2),:));
 
   [CC]=corrcoef(A1,A2,'rows','complete');
    allCmaxW(ind1,ind2)= CC(1,2);%max(abs(c));
   ;
end;end
    allCmaxW(allCmaxW>0.99)=NaN;
      allCmaxW(allCmaxW==0)=NaN;
    
    cindi= [cindi; allCmax(:)];
cwide= [cwide; allCmaxW(:)];
rdist= [rdist; ROIdist(:)];

    %%
     subthresIndi= indifile.allresults.roaster;
   subthresWide= widefile.allresults.roaster;
rate_thres=5;
    clear allCmax  ROIdist
for ind1=1:length(nsel)
for ind2=1:length(nsel)
  A1=  zscore(subthresIndi(nsel(ind1),:));A2=  zscore(subthresIndi(nsel(ind2),:));
%  if ind1<ind2
 if  length(find(A1>0)) >rate_thres & length(find(A2>0)) >rate_thres
    [c,lags]= xcorr(fastsmooth(A1,10,1,1),fastsmooth(A2,10,1,1),100,'Coeff');
    c1=c(95:105);
    mean_corr = mean(c1);
    allCmax(ind1,ind2)= mean_corr;
    allCx= [allCx ; c];
    allCm=[allCm; mean_corr];
 %     allCx= [allCx ; c];
 %allCmax(ind2,ind1)= c1(n2);

 else;  allCmax(ind1,ind2)=NaN;end;
   ROIdist(ind1,ind2)= sqrt(sum((irm(:,nsel(ind1)).*sc-irm(:,nsel(ind2)).*sc).^2));
   % ROIdist(ind2,ind1)= ROIdist(ind1,ind2);m
%  end
end;end
    allCmax(allCmax>0.99)=NaN;
      allCmax(allCmax==0)=NaN;
    
  clear allCmaxW  
for ind1=1:length(nsel)
for ind2=1:length(nsel)
  A1=  zscore(subthresWide(nsel(ind1),:));A2=  zscore(subthresWide(nsel(ind2),:));
    if ind1<ind2
   if  length(find(A1>0)) >rate_thres & length(find(A2>0)) >rate_thres
    [c,lags]= xcorr(fastsmooth(A1,10,1,1),fastsmooth(A2,10,1,1),100,'Coeff');
    c1=c(95:105); % This is the +/- around 0 from the 200 parameter, the offset
    mean_corr = mean(c1); % Use mean
    allCmaxW(ind1,ind2)= mean_corr;
     allCmaxW(ind2,ind1)= mean_corr;
   % allCmaxW(ind2,ind1)=  allCmax(ind1,ind2);
   else;  allCmaxW(ind1,ind2)=NaN;end
    end
end;end
    allCmaxW(allCmaxW>0.99)=NaN;  
         allCmaxW(allCmaxW==0)=NaN;
    
    
cindiS= [cindiS; allCmax(:)];
cwideS= [cwideS; allCmaxW(:)];


  
end
rdist2=rdist;
cindi=cindi(~isnan(cwide));
rdist=rdist(~isnan(cwide));
cwide=cwide(~isnan(cwide));

cindiS=cindiS(~isnan(cwideS));
%allCx=allCx(~isnan(cwideS),:);
rdist2=rdist2(~isnan(cwideS));
cwideS=cwideS(~isnan(cwideS));
rdist2=rdist2(~isnan(cindiS));
cwideS=cwideS(~isnan(cindiS));
cindiS=cindiS(~isnan(cindiS));
%allCx=allCx(~isnan(cindiS),:)

% Ignore very close neurons. This was a condition for selecting matching neurons between indi and wide 
use_idx = find(rdist2 > 8); % Greater than 8um between neuron pairs
rdist2 = rdist2(use_idx);
cindiS = cindiS(use_idx);
cwideS = cwideS(use_idx);

use_idx = find(rdist > 8);
rdist = rdist(use_idx);
cindi = cindi(use_idx);
cwide = cwide(use_idx);


%% Plot all of the data
% PLOT subthreshold cross correlation with regression line
figure('Color','w', 'Renderer', 'painters', 'Position', [0 0 900 800]);
plot(rdist,cindi,'.r','Markersize',10)
hold on, plot(rdist,cwide,'.b','Markersize',10)
fitResults1 = polyfit(rdist,cindi,1);
yplot1 = polyval(fitResults1,rdist);
plot(rdist,yplot1,'-r','Linewidth',2)
fitResults2 = polyfit(rdist,cwide,1);
yplot1 = polyval(fitResults2,rdist);
plot(rdist,yplot1,'-b','Linewidth',2)
ylabel('Corr')
legend({['indi slope=' num2str(fitResults1(1))], ['wide slope=' num2str(fitResults2(1))]});

title_string = ['Cross correlation subthreshold (Vm) indi vs. wide'];
title(title_string);
saveas(gcf, [save_fig_path 'Cross Correlation/Jpeg Format/' title_string '.jpg']);
saveas(gcf, [save_fig_path 'Cross Correlation/SVG Format/' title_string '.svg']);
saveas(gcf, [save_fig_path 'Cross Correlation/EPS Format/' title_string '.eps'], 'epsc');


% Plot subthreshold cross correlation over binned distance
bin_size = 30; % In um
last_bin = 180; % Minimum value of the last bin so that there are enough points
indi_corr_bins = [];
wide_corr_bins = [];
ft_corr_cells = {};

all_bins = [];
bin_labels = {};
for i=0:bin_size:last_bin
    
    % Bin all pairs that are greater than 180 um
    if i >= last_bin
        corr_idx = find(rdist >= i);
    else
        corr_idx = find(rdist >= i & rdist < i+bin_size);
    end
    
    indi_corr_bins = horzcat_pad(indi_corr_bins, cindi(corr_idx));
    wide_corr_bins = horzcat_pad(wide_corr_bins, cwide(corr_idx));
    
    % Store the correlations in bins using cell arrays
    ft_corr_cells = cat(1, [num2cell(cindi(corr_idx)), num2cell(cwide(corr_idx))]);
    
    % Store all of the correlation distributions in one matrix
    all_bins = horzcat_pad(all_bins, cindi(corr_idx));
    all_bins = horzcat_pad(all_bins, cwide(corr_idx));
    
    if i >= last_bin
        bin_labels = cat(2, bin_labels, ['>=' num2str(i)]);
    else
        bin_labels = cat(2, bin_labels, [num2str(i) '-' num2str(i+bin_size - 1)]);
    end

    bin_labels = cat(2, bin_labels, [' ']);
end
sub_indi_corr_bin = indi_corr_bins;
sub_wide_corr_bin = wide_corr_bins;

% X = categorical(bin_labels);
% X = reordercats(X, bin_labels);
figure('Renderer', 'painters', 'Position', [0 0 900 800]);
boxplot(all_bins, 'labels', bin_labels, 'notch', 'on', 'colors', [ 0.4 0.4 0.4], 'symbol', '.k');
%boxplotGroup({indi_corr_bins, wide_corr_bins}, 'PrimaryLabels', {'indi', 'wide'},...
%        'SecondaryLabels', bin_labels);

%a = get(get(gca,'children'),'children');
%t = get(a,'tag')

ylabel('Pearson''s correlation');
title_string = ['Subthreshold Vm cross correlation over binned distance ' num2str(bin_size)];
title(title_string);

saveas(gcf, [save_fig_path 'Cross Correlation/Jpeg Format/' title_string '.jpg']);
saveas(gcf, [save_fig_path 'Cross Correlation/SVG Format/' title_string '.svg']);
saveas(gcf, [save_fig_path 'Cross Correlation/EPS Format/' title_string '.eps'], 'epsc');

% Statistical test for subthreshold correlations from binned distances
% Will keep the columns as DMD and wide field, because that is the factor
% we are looking for interaction between the individual and wide
% I think I have to line up the same bins across both conditions, add Nans
% intermittently. This may allow me to do the Friedman's maybe?
ft_corr = [];
max_reps = 0;
for i=1:size(indi_corr_bins, 2)
    ft_corr = [ft_corr; horzcat_pad(nanmean(indi_corr_bins(:, i)), nanmean(wide_corr_bins(:, i)))];
end

%-- Generate stats table for Vm-Vm binned distance correlations --
vm_stats = [];
vm_bin_stats = ["Binned Distance", "Median Targeted", "Median Widefield", "Number of neuron Pairs"];
for i=1:(length(bin_labels)/2)
    vm_bin_stats = [vm_bin_stats; ...
        bin_labels(2*i - 1), nanmedian(all_bins(:, 2*i - 1)), nanmedian(all_bins(:, 2*i)), sum(~isnan(all_bins(:, 2*i - 1)))];
end
vm_stats = [vm_stats, vm_bin_stats];

%  Kolmogorov-Smirnov test for difference in the two correlations over
%  distance
disp('Subthreshold correlation over distance (Comparing DMD vs. Wide Field) Kolmogorov–Smirnov test:');
[ordered_dis, sidx ] = sort(rdist);
sorted_cindi = cindi(sidx);
sorted_cwide = cwide(sidx);

[h, p, ks2stat] = kstest2(sorted_cindi, sorted_cwide)

% Add kstest output to vm_stats_table
vm_stats = [vm_stats, ["Kolmogorov-Smirnov test", ""; "P-value", "Statistic"; p, ks2stat; repmat("", 5, 2)]];

% Unfortunantely it does not like the cell arrays
%disp('Friedman''s test on Subthreshold correlations:');
%[p,tbl,stats] = friedman(ft_corr_cells, 1)

disp('Friedman''s test on Subthreshold correlations DMD in column 1, wide field in column 2 (using averages):');
[p,tbl,stats] = friedman(ft_corr, 1)

vm_stats = [vm_stats, ...
    ["Friedman's test", repmat("", 1, size(tbl, 2) - 1)   ; ...
    tbl]]; 

disp('Kruskal-Wallis test on subthreshold DMD:');
[p,tbl,stats] = kruskalwallis(indi_corr_bins)

vm_stats = [vm_stats, ...
    ["Kruskal-Wallis for targeted", repmat("", 1, size(tbl, 2) - 1)   ; ...
    tbl]]; 

disp('Kruskal-Wallis test on subthreshold wide field:');
[p,tbl,stats] = kruskalwallis(wide_corr_bins)

vm_stats = [vm_stats, ...
    ["Kruskal-Wallis for Widefield", repmat("", 1, size(tbl, 2) - 1)   ; ...
    tbl]]; 

writematrix(vm_stats, [save_fig_path 'Cross Correlation/Data Tables/vm_binned_stats.csv']);

% Plot spike-spike correlations over binned distances
indi_corr_bins = [];
wide_corr_bins = [];
all_bins = [];
bin_labels = {};
for i=0:bin_size:last_bin
    if i >= last_bin
        corr_idx = find(rdist2 >= i);
    else
        corr_idx = find(rdist2 >= i & rdist2 < i+bin_size);
    end
    indi_corr_bins = horzcat_pad(indi_corr_bins, cindiS(corr_idx));
    wide_corr_bins = horzcat_pad(wide_corr_bins, cwideS(corr_idx));
    
    % Store each FOVs distribution as column vectors
    all_bins = horzcat_pad(all_bins, cindiS(corr_idx));
    all_bins = horzcat_pad(all_bins, cwideS(corr_idx));
    
    
    if i >= last_bin
        bin_labels = cat(2, bin_labels, ['>=' num2str(i)]);
    else
        bin_labels = cat(2, bin_labels, [num2str(i) '-' num2str(i+bin_size - 1)]);
    end

    bin_labels = cat(2, bin_labels, ['-']);
    
end

ft_corr = [];
max_reps = 0;
for i=1:size(indi_corr_bins, 2)
    ft_corr = [ft_corr; horzcat_pad(nanmean(indi_corr_bins(:, i)), nanmean(wide_corr_bins(:, i)))];
end

% Spike-spike statistical tests
disp('Friedman''s test on spike-spike correlations DMD in column 1, wide field in column 2 (using averages):');
[p,tbl,stats] = friedman(ft_corr, 1)

disp('Kruskal-Wallis test on DMD spike-spike correlations over binned distance:');
[p,tbl,stats] = kruskalwallis(indi_corr_bins)

disp('Kruskal-Wallis test on wide field spike-spike correlations over binned distance:');
[p,tbl,stats] = kruskalwallis(wide_corr_bins)

%X = categorical(bin_labels);
%X = reordercats(X, bin_labels);

figure('Renderer', 'painters', 'Position', [0 0 900 800]);
boxplot(all_bins, 'labels', bin_labels, 'notch', 'on', 'colors', [ 0.4 0.4 0.4], 'symbol', '.k');
%boxplotGroup({indi_corr_bins, wide_corr_bins}, 'PrimaryLabels', {'indi', 'wide'}, 'SecondaryLabels', bin_labels);
ylabel('Pearson''s correlation');
title_string = ['Spike-Spike cross correlation over binned distance ' num2str(bin_size)];
title(title_string);

saveas(gcf, [save_fig_path 'Cross Correlation/Jpeg Format/' title_string '.jpg']);
saveas(gcf, [save_fig_path 'Cross Correlation/EPS Format/' title_string '.eps'], 'epsc');
saveas(gcf, [save_fig_path 'Cross Correlation/SVG Format/' title_string '.svg']);



%-- Generate stats table for spike-spike binned distance correlations --
% Spike-Spike correlation bins
ss_stats_table = ["Binned Distance", "Median Targeted", "Median Widefield", "Number of neuron Pairs"];
for i=1:(length(bin_labels)/2)
    ss_stats_table = [ss_stats_table; ...
        bin_labels(2*i - 1), nanmedian(all_bins(:, 2*i - 1)), nanmedian(all_bins(:, 2*i)), sum(~isnan(all_bins(:, 2*i - 1)))];
end

writematrix(ss_stats_table, [save_fig_path 'Cross Correlation/Data Tables/spike_binned_stats.csv']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spike to spike correlation

figure('Color','w', 'Renderer', 'painters', 'Position', [0 0 900 800])
plot(rdist2,cindiS,'.r','Markersize',10)
hold on,plot(rdist2,cwideS,'.b','Markersize',10)
fitResults1 = polyfit(rdist2,cindiS,1);
yplot1 = polyval(fitResults1,rdist2);
plot(rdist2,yplot1,'-r','Linewidth',2)
fitResults2 = polyfit(rdist2,cwideS,1);
yplot1 = polyval(fitResults2,rdist2);
plot(rdist2,yplot1,'-b','Linewidth',2)
xlabel('Distance (\mum)');
ylim([-0.05, 0.55]);
ylabel('Corr');
legend({['indi slope=' num2str(fitResults1(1))], ['wide slope=' num2str(fitResults2(1))]});
title_string = ['Cross correlation spike-spike indi vs. wide'];
title(title_string);

saveas(gcf, [save_fig_path 'Cross Correlation/Jpeg Format/' title_string '.jpg']);
saveas(gcf, [save_fig_path 'Cross Correlation/EPS Format/' title_string '.eps'], 'epsc');
saveas(gcf, [save_fig_path 'Cross Correlation/SVG Format/' title_string '.svg']);

%  Kolmogorov-Smirnov test for difference in the two correlations over
%  distance
disp('Spike-spike correlation over distance KS test:');
[ordered_dis, sidx ] = sort(rdist2);
sorted_cindiS = cindiS(sidx);
sorted_cwideS = cwideS(sidx);

[h, p, ks2stat] = kstest2(sorted_cindiS, sorted_cwideS)

%---- Compare the Vm-Vm correlation with the spike-spike correlations
disp('Vm vs spike correlations for widefield:');
[h, p, ci, stats] = ttest2(cwide, cwideS)


disp('Vm vs spike correlations for targeted:');
[h, p, ci, stats] = ttest2(cindi, cindiS)

