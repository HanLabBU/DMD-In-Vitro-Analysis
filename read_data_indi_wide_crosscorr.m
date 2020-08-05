
% In vivo path
%cd('\\engnas.bu.edu\research\eng_research_handata\EricLowet\DMD\invivoDMD\')
%cd('\\engnas.bu.edu\research\eng_research_handata\EricLowet\DMD\comp_wide_indi\')

close all;
clear all;

addpath('.');

% In vitro path
cd('\\engnas.bu.edu\research\eng_research_handata\Pierre Fabris\DMD Project\All In Vitro Analysis\');

% Final figure folder
save_fig_path = '\\engnas.bu.edu\research\eng_research_hanlab\DMD Paper\In Vitro Plots\';

% Scripts for needed functions
addpath('\\ad\eng\research\eng_research_handata\EricLowet\DMD\main_analysis\');

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

 sc=0.1625*2; % micrometer per pixel 40x
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
   % ROIdist(ind2,ind1)= ROIdist(ind1,ind2);
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

cindiS=cindiS(~isnan(cwideS))
%allCx=allCx(~isnan(cwideS),:)
rdist2=rdist2(~isnan(cwideS))
cwideS=cwideS(~isnan(cwideS))
rdist2=rdist2(~isnan(cindiS))
cwideS=cwideS(~isnan(cindiS))
cindiS=cindiS(~isnan(cindiS))
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
% PLOT subthreshold cross correlation
figure('Color','w', 'Renderer', 'painters')
plot(rdist,cindi,'.r','Markersize',10)
hold on, plot(rdist,cwide,'.k','Markersize',10)
fitResults1 = polyfit(rdist,cindi,1);
yplot1 = polyval(fitResults1,rdist);
plot(rdist,yplot1,'-r','Linewidth',2)
fitResults2 = polyfit(rdist,cwide,1);
yplot1 = polyval(fitResults2,rdist);
plot(rdist,yplot1,'-k','Linewidth',2)
ylabel('Corr')
legend({['indi slope=' num2str(fitResults1(1))], ['wide slope=' num2str(fitResults2(1))]});

title_string = ['Cross correlation subthreshold (Vm) indi vs. wide'];
title(title_string);
saveas(gcf, [save_fig_path 'Cross Correlation/Jpeg Format/' title_string '.jpg']);
saveas(gcf, [save_fig_path 'Cross Correlation/SVG Format/' title_string '.svg']);
saveas(gcf, [save_fig_path 'Cross Correlation/EPS Format/' title_string '.eps'], 'epsc');

% Plot subthreshold cross correlation over binned distance
bin_size = 50; % In um
last_bin = 180;
indi_corr_bins = [];
wide_corr_bins = [];
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
    
    % Store all of the correlation distributions in one matrix
    all_bins = horzcat_pad(all_bins, cindi(corr_idx));
    all_bins = horzcat_pad(all_bins, cwide(corr_idx));
    
    if i >= last_bin
        bin_labels = cat(2, bin_labels, ['>=' num2str(i)]);
    else
        bin_labels = cat(2, bin_labels, [num2str(i) '-' num2str(i+bin_size - 1)]);
    end
end
X = categorical(bin_labels);
X = reordercats(X, bin_labels);
figure('Renderer', 'painters', 'Position', [0 0 900 800]);
boxplotGroup({indi_corr_bins, wide_corr_bins}, 'PrimaryLabels', {'indi', 'wide'},...
        'SecondaryLabels', bin_labels);
ylabel('Pearson''s correlation');
title_string = ['Subthreshold Vm cross correlation over binned distance ' num2str(bin_size)];
title(title_string);

saveas(gcf, [save_fig_path 'Cross Correlation/Jpeg Format/' title_string '.jpg']);
saveas(gcf, [save_fig_path 'Cross Correlation/SVG Format/' title_string '.svg']);
saveas(gcf, [save_fig_path 'Cross Correlation/EPS Format/' title_string '.eps'], 'epsc');

%--Perform shuffling of the data to test significance in the calculated regression lines
num_reshuffles = 1000;
indi_shuf_slopes = [];
wide_shuf_slopes = [];
slope_diff = [];

% Shuffle subthreshold trace identities
dist_corr_pairs = [cindi'; cwide']; 
for i=1:num_reshuffles
    % Assign random individual correlations
    % Logical indexing is by column, so the cindis are transposed
    indi_idx = logical(randi([0 1], [1 length(rdist)]));
    indi_idx = [indi_idx; ~indi_idx];
    wide_idx = ~indi_idx;
       
	indi_shuff_cor = dist_corr_pairs(indi_idx);
    wide_shuff_cor = dist_corr_pairs(wide_idx);
        
    % Calculate best fit slope and store
    indi_fit = polyfit(rdist, indi_shuff_cor,1);
    wide_fit = polyfit(rdist, wide_shuff_cor,1);
    
    indi_shuf_slopes = [indi_shuf_slopes, indi_fit(1)];
    wide_shuf_slopes = [wide_shuf_slopes, wide_fit(1)];
    
    % Store the slope difference indi - wide
    slope_diff = [slope_diff, indi_fit(1) - wide_fit(1)];
end

indi_top = nanmean(indi_shuf_slopes) + 3*nanstd(indi_shuf_slopes)
indi_bot = nanmean(indi_shuf_slopes) - 3*nanstd(indi_shuf_slopes)

wide_top = nanmean(wide_shuf_slopes) + 3*nanstd(wide_shuf_slopes)
wide_bot = nanmean(wide_shuf_slopes) - 3*nanstd(wide_shuf_slopes)

% Plot the shuffled distributions of the correlations over distance
figure('Renderer', 'painters', 'Position', [0 0 1700 700]);
subplot(1, 2, 1);
xline(fitResults1(1), 'LineWidth', 2, 'Color', 'red');
hold on;
xline(nanmean(indi_shuf_slopes), 'LineWidth', 2, 'Color', 'black');
hold on;
xline(indi_top, 'LineWidth', 2, 'Color', 'green');
hold on;
xline(indi_bot, 'LineWidth', 2, 'Color', 'green');
[counts, bin_centers] = hist(indi_shuf_slopes);
bar(bin_centers, counts, 'BarWidth', 1);
legend({'Observed','Mean shuffled', '\mu+/- 3sigma'});
title('Individual Shuffled correlations');

subplot(1, 2, 2);
[counts, bin_centers] = hist(wide_shuf_slopes);
xline(fitResults2(1), 'LineWidth', 2, 'Color', 'red');
hold on;
xline(nanmean(wide_shuf_slopes), 'LineWidth', 2, 'Color', 'black');
hold on;
xline(wide_top, 'LineWidth', 2, 'Color', 'green');
hold on;
xline(wide_bot, 'LineWidth', 2, 'Color', 'green');
hold on;
bar(bin_centers, counts, 'BarWidth', 1);
legend({'Observed','Mean shuffled', '\mu+/- 3sigma'});
title('Wide Field Shuffled correlations');

sgtitle('Shuffling subthreshold Vm correlations');

% Plot the shuffled slope differences
diff_top = nanmean(slope_diff) + 3*nanstd(slope_diff);
diff_bot = nanmean(slope_diff) - 3*nanstd(slope_diff);
figure('Renderer', 'painters');
xline(fitResults1(1) - fitResults2(1), 'LineWidth', 2, 'Color', 'red');
hold on;
xline(nanmean(slope_diff), 'LineWidth', 2, 'Color', 'black');
hold on;
xline(diff_top, 'LineWidth', 2, 'Color', 'green');
hold on;
xline(diff_bot, 'LineWidth', 2, 'Color', 'green');
hold on;
[counts, bin_centers] = hist(slope_diff);
bar(bin_centers, counts, 'BarWidth', 1);
legend({'Observed','Mean shuffled', '\mu+/- 3sigma'});
title('Shuffled Slope Difference From Subthreshold Vm (Individual - Wide Field)');

%TODO save the shuffled plots

%{
% Plot subthreshold cross correlations
[h,p,ci,stats] = ttest(cindi,cwide)
figure('COlor','w','Position', [ 300 300 200 200])
V1=nanmean(cindi);V1s=(std(cindi)./sqrt(length(cindi)));
V2=nanmean(cwide);V2s=(std(cindi)./sqrt(length(cwide)));
bar( [ 1 ], [V1],0.7,'FaceColor', [ 0.7 0.2 0.1]) , hold on,bar( [ 2 ], [V2],0.7,'FaceColor', [0.1 0.4 0.7])
set(gca,'Xtick', [ 1 2],'Xticklabel', {'DMD' ; 'Widefield'})
errorbar([ 1 2], [ V1 V2], [V1s V2s],'.k','Linewidth', 2)
axis tight;ylabel('corr max')
xlim([ 0.5 2.5]); %ylim([0  20])
title([ 'p= ' num2str(p)])
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spike to spike correlation

figure('Color','w', 'Renderer', 'painters')
plot(rdist2,cindiS,'.r','Markersize',10)
hold on,plot(rdist2,cwideS,'.k','Markersize',10)
fitResults1 = polyfit(rdist2,cindiS,1);
yplot1 = polyval(fitResults1,rdist2);
plot(rdist2,yplot1,'-r','Linewidth',2)
fitResults2 = polyfit(rdist2,cwideS,1);
yplot1 = polyval(fitResults2,rdist2);
plot(rdist2,yplot1,'-k','Linewidth',2)
xlabel('Distance (\mum)');
ylabel('Corr');
legend({['indi slope=' num2str(fitResults1(1))], ['wide slope=' num2str(fitResults2(1))]});
title_string = ['Cross correlation spike-spike indi vs. wide'];
title(title_string);

saveas(gcf, [save_fig_path 'Cross Correlation/Jpeg Format/' title_string '.jpg']);
saveas(gcf, [save_fig_path 'Cross Correlation/EPS Format/' title_string '.eps'], 'epsc');
saveas(gcf, [save_fig_path 'Cross Correlation/SVG Format/' title_string '.svg']);

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
    
end
X = categorical(bin_labels);
X = reordercats(X, bin_labels);

figure('Renderer', 'painters', 'Position', [0 0 900 800]);
boxplotGroup({indi_corr_bins, wide_corr_bins}, 'PrimaryLabels', {'indi', 'wide'}, 'SecondaryLabels', bin_labels);
ylabel('Pearson''s correlation');
title_string = ['Spike-Spike cross correlation over binned distance ' num2str(bin_size)];
title(title_string);

saveas(gcf, [save_fig_path 'Cross Correlation/Jpeg Format/' title_string '.jpg']);
saveas(gcf, [save_fig_path 'Cross Correlation/EPS Format/' title_string '.eps'], 'epsc');
saveas(gcf, [save_fig_path 'Cross Correlation/SVG Format/' title_string '.svg']);

%--Perform shuffling of the data to test significance in the calculated regression lines
indi_shuf_slopes = [];
wide_shuf_slopes = [];
slope_diff = [];

% Shuffle spike-spike correlation identities
dist_corr_pairs = [cindiS'; cwideS']; 
for i=1:num_reshuffles
    % Assign random individual correlations
    % Logical indexing is by column, so the cindis are transposed
    indi_idx = logical(randi([0 1], [1 length(rdist2)]));
    indi_idx = [indi_idx; ~indi_idx];
    wide_idx = ~indi_idx;
       
	indi_shuff_cor = dist_corr_pairs(indi_idx);
    wide_shuff_cor = dist_corr_pairs(wide_idx);
        
    % Calculate best fit slope and store
    indi_fit = polyfit(rdist2, indi_shuff_cor,1);
    wide_fit = polyfit(rdist2, wide_shuff_cor,1);
    
    indi_shuf_slopes = [indi_shuf_slopes, indi_fit(1)];
    wide_shuf_slopes = [wide_shuf_slopes, wide_fit(1)];
    
    % Store the slope difference indi - wide
    slope_diff = [slope_diff, indi_fit(1) - wide_fit(1)];
end

indi_top = nanmean(indi_shuf_slopes) + 3*nanstd(indi_shuf_slopes)
indi_bot = nanmean(indi_shuf_slopes) - 3*nanstd(indi_shuf_slopes)

wide_top = nanmean(wide_shuf_slopes) + 3*nanstd(wide_shuf_slopes)
wide_bot = nanmean(wide_shuf_slopes) - 3*nanstd(wide_shuf_slopes)

% Plot the shuffled distributions of the correlations over distance
figure('Renderer', 'painters', 'Position', [0 0 1700 700]);
subplot(1, 2, 1);
xline(fitResults1(1), 'LineWidth', 2, 'Color', 'red');
hold on;
xline(nanmean(indi_shuf_slopes), 'LineWidth', 2, 'Color', 'black');
hold on;
xline(indi_top, 'LineWidth', 2, 'Color', 'green');
hold on;
xline(indi_bot, 'LineWidth', 2, 'Color', 'green');
[counts, bin_centers] = hist(indi_shuf_slopes);
bar(bin_centers, counts, 'BarWidth', 1);
legend({'Observed','Mean shuffled', '\mu+/- 3sigma'});
title('Individual Shuffled correlations');

subplot(1, 2, 2);
[counts, bin_centers] = hist(wide_shuf_slopes);
xline(fitResults2(1), 'LineWidth', 2, 'Color', 'red');
hold on;
xline(nanmean(wide_shuf_slopes), 'LineWidth', 2, 'Color', 'black');
hold on;
xline(wide_top, 'LineWidth', 2, 'Color', 'green');
hold on;
xline(wide_bot, 'LineWidth', 2, 'Color', 'green');
hold on;
bar(bin_centers, counts, 'BarWidth', 1);
legend({'Observed','Mean shuffled', '\mu+/- 3sigma'});
title('Wide Field Shuffled correlations');

sgtitle('Shuffling Spike-Spike correlations');

% Plot the shuffled slope differences
diff_top = nanmean(slope_diff) + 3*nanstd(slope_diff);
diff_bot = nanmean(slope_diff) - 3*nanstd(slope_diff);
figure('Renderer', 'painters');
xline(fitResults1(1) - fitResults2(1), 'LineWidth', 2, 'Color', 'red');
hold on;
xline(nanmean(slope_diff), 'LineWidth', 2, 'Color', 'black');
hold on;
xline(diff_top, 'LineWidth', 2, 'Color', 'green');
hold on;
xline(diff_bot, 'LineWidth', 2, 'Color', 'green');
hold on;
[counts, bin_centers] = hist(slope_diff);
bar(bin_centers, counts, 'BarWidth', 1);
legend({'Observed','Mean shuffled', '\mu+/- 3sigma'});
title('Shuffled Slope Difference from spike-spike (Individual - Wide Field)');

%TODO save the shuffled plots


%{
[h,p,ci,stats] = ttest(cindiS,cwideS)
figure('COlor','w','Position', [ 300 300 200 200])
V1=nanmean(cindiS);V1s=(std(cindiS)./sqrt(length(cindiS)));
V2=nanmean(cwideS);V2s=(std(cindiS)./sqrt(length(cwideS)));
bar( [ 1 ], [V1],0.7,'FaceColor', [ 0.7 0.2 0.1]) , hold on,bar( [ 2 ], [V2],0.7,'FaceColor', [0.1 0.4 0.7])
set(gca,'Xtick', [ 1 2],'Xticklabel', {'DMD' ; 'Widefield'})
errorbar([ 1 2], [ V1 V2], [V1s V2s],'.k','Linewidth', 2)
axis tight;ylabel('corr max')
xlim([ 0.5 2.5]); %ylim([0  20])
title([ 'Spike-Spike bar graph p= ' num2str(p)])
%}



% %% PLOT
% figure('COlor','w'),plot(indiB,'r'); hold on,plot(wideB,'k')
% legend indi wide
% [h,p,ci,stats] = ttest(indiB,wideB)
% figure('COlor','w','Position', [ 300 300 200 200])
% V1=((1-nanmean(indiB)).*-1).*100;V1s=(std(indiB)./sqrt(length(indiB))).*100;
% V2=((1-nanmean(wideB)).*-1).*100;V2s=(std(wideB)./sqrt(length(wideB))).*100;
% bar( [ 1 ], [V1],0.7,'FaceColor', [ 0.7 0.2 0.1]) , hold on,bar( [ 2 ], [V2],0.7,'FaceColor', [0.1 0.4 0.7])
% set(gca,'Xtick', [ 1 2],'Xticklabel', {'DMD' ; 'Widefield'})
% errorbar([ 1 2], [ V1 V2], [V1s V2s],'.k','Linewidth', 2)
% axis tight;ylabel('signal reduction %')
% xlim([ 0.5 2.5]); %ylim([0  20])
% title([ 'p= ' num2str(p)])



