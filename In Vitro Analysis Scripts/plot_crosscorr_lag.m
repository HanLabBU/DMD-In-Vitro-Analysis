clc;
clear all;
close all;

% Maingear office computer directory
cd('~/handata_server/Pierre Fabris/DMD Project/All In Vitro Analysis');
%cd('\\engnas.bu.edu\research\eng_research_handata\EricLowet\DMD\comp_wide_indi\')

save_fig_path = '~/handata_server/Pierre Fabris/DMD Project/Data Figures/';

addpath(genpath('~/handata_server/EricLowet/Scripts/'));

ses=dir('*.mat');

%%% search for indi
clear findwide
for id=1:length(ses)
    if   strfind(ses(id).name, 'Wide')>0
        findwide(id)=1;
    else
        findwide(id)=0;
    end
    
end

wideloc=find(findwide)
indiloc=find([~findwide]);

cindi=[];cwide=[]; rdist=[];allR=[];
cindiS=[];cwideS=[]; allCx=[];allCm=[];
rdist2=[];allRW=[]; 
allCxW=[];allCmW=[];


%Calculate cross correlation for individual DMD
for id=1:length(indiloc)
	F=  ses(indiloc(id)).name;
    	%   widefile=load(ses(wideloc(id)).name);
    	indifile=load(ses(indiloc(id)).name);
    	%     clear wrm % ROI centroid
    	%     for id2= 1:length(widefile.allresults.roi)
    	%         [ x y]=find(widefile.allresults.roi{id2});
    	%         wrm(:,id2)= round(mean([x , y]));end
    	clear irm
    	for id2= 1:length(indifile.allresults.roi)
        	[ x y]=find(indifile.allresults.roi{id2});
        	irm(:,id2)= round(mean([x , y]));
	end
 

	nsel=1:size(irm,2);sc=0.1625*2; % number of microns per pixel

	%mROI=[];
    	%for id3=1:size(irm,2) % matching ROI
    	%    cents=irm(:,id3);
    	%    pxdiff=(sqrt(sum(bsxfun(@minus, wrm, cents).^2)));
    	%
    	%  wloc=find(pxdiff<8);
    	%  if ~isempty(wloc)
    	%mROI=[mROI; [ id3 wloc]];end

	subthresIndi= indifile.allresults.orig_traceDN;
 	%subthresWide= widefile.allresults.orig_traceDN;
	
	%for ix=1:size( subthresIndi,1)
	%  subthresIndi(ix,:)= subthresIndi(ix,:)-fastsmooth( subthresIndi(ix,:),400,1,1);
	%end
   
	%subthresWide= subthresWide-fastsmooth( subthresWide,10,1,1);
	
	clear allCmax  ROIdist

	for ind1=1:length(nsel)
		for ind2=1:length(nsel)
   			if ind1<ind2
				% A1=  zscore(subthresIndi(nsel(ind1),:));A2=  zscore(subthresIndi(nsel(ind2),:));
   				A1=  (subthresIndi(nsel(ind1),:));A2=  (subthresIndi(nsel(ind2),:));
   				[CC]=corrcoef(A1,A2,'Rows','complete');
 				% [c,lags]= xcorr(A1,A2,10,'Coeff');
 		  		allCmax(ind1,ind2)= CC(1,2);%max(abs(c));
 		   		allCmax(ind2,ind1)=allCmax(ind1,ind2);end
 		 		ROIdist(ind1,ind2)= sqrt(sum((irm(:,nsel(ind1)).*sc -irm(:,nsel(ind2)).*sc).^2));
 		 		%ROIdist(ind1,ind2)=ROIdist(ind2,ind1);
			end
		end
   
	allCmax(allCmax>0.99)=NaN;
	allCmax(allCmax==0)=NaN;
   
	cindi= [cindi; allCmax(:)];
	%cwide= [cwide; allCmaxW(:)];
	rdist= [rdist; ROIdist(:)];

	%% 
	
	subthresIndi= indifile.allresults.roaster;
	%subthresWide= widefile.allresults.roaster;

	clear allCmax  ROIdist
	for ind1=1:length(nsel)
		for ind2=1:length(nsel)
	   		%if ind1<ind2
	  		A1=  zscore(subthresIndi(nsel(ind1),:));A2=  zscore(subthresIndi(nsel(ind2),:));
	 		
			if  length(find(A1>0)) >5 & length(find(A2>0)) >5
				%[c,lags]= xcorr(fastsmooth(A1,10,1,1),fastsmooth(A2,10,1,1),10,'Coeff');
				%allCmax(ind1,ind2)= max(abs(c));
	         		[c,lags]= xcorr(fastsmooth(A1,10,1,1),fastsmooth(A2,10,1,1),200,'Coeff');
	    			c1=c(191:211);
	    			[n1 n2]= max(abs(c1));
	     			
				if ind1~=ind2
	    				allCmax(ind1,ind2)= mean(c1);
	    				allCx= [allCx ; c];
	    				allCm=[allCm;c1(n2)];
	    				allR= [allR; sqrt(sum((irm(:,nsel(ind1))-irm(:,nsel(ind2))).^2))];
	     			else
	         			allCmax(ind1,ind2)= NaN;
	          			allCx= [allCx ; c.*NaN];
	            			allCm=[allCm;NaN];
	                		allR= [allR; NaN];
	   			end
	    
	 		else
				allCmax(ind1,ind2)=NaN;
			end
	  		
			% ROIdist(ind1,ind2)= sqrt(sum((irm(:,nsel(ind1))-irm(:,nsel(ind2))).^2));
		end
	end

	allCmax(allCmax>0.99)=NaN; %0.35 V1

	%   clear allCmaxW  
	% for ind1=1:length(nsel)
	% for ind2=1:length(nsel)
	%   A1=  zscore(subthresWide(nsel(ind1),:));A2=  zscore(subthresWide(nsel(ind2),:));
	%    if  length(find(A1>0)) >5 & length(find(A2>0)) >5
	%      [c,lags]= xcorr(fastsmooth(A1,10,1,1),fastsmooth(A2,10,1,1),10,'Coeff');
	%     allCmaxW(ind1,ind2)= max(abs(c)); else;  allCmaxW(ind1,ind2)=NaN;end
	% 
	% end;end
	%     allCmaxW(allCmaxW>0.6)=NaN;  
	%     
	    
	    
	cindiS= [cindiS; allCmax(:)];
	%cwideS= [cwideS; allCmaxW(:)];
  
end


%% Perform the lag cross correlation for the wide field only
for id=1:length(wideloc)
	F=  ses(wideloc(id)).name;
    	%   widefile=load(ses(wideloc(id)).name);
    	widefile=load(ses(wideloc(id)).name);
    	%     clear wrm % ROI centroid
    	%     for id2= 1:length(widefile.allresults.roi)
    	%         [ x y]=find(widefile.allresults.roi{id2});
    	%         wrm(:,id2)= round(mean([x , y]));end
    	clear irm
    	for id2= 1:length(widefile.allresults.roi)
        	[ x y]=find(widefile.allresults.roi{id2});
        	irm(:,id2)= round(mean([x , y]));
	end
 

	nsel=1:size(irm,2);sc=0.1625*2; % number of microns per pixel

	%mROI=[];
    	%for id3=1:size(irm,2) % matching ROI
    	%    cents=irm(:,id3);
    	%    pxdiff=(sqrt(sum(bsxfun(@minus, wrm, cents).^2)));
    	%
    	%  wloc=find(pxdiff<8);
    	%  if ~isempty(wloc)
    	%mROI=[mROI; [ id3 wloc]];end

	subthresIndi= widefile.allresults.orig_traceDN;
 	%subthresWide= widefile.allresults.orig_traceDN;
	
	%for ix=1:size( subthresIndi,1)
	%  subthresIndi(ix,:)= subthresIndi(ix,:)-fastsmooth( subthresIndi(ix,:),400,1,1);
	%end
   
	%subthresWide= subthresWide-fastsmooth( subthresWide,10,1,1);
	
	clear allCmaxW  ROIdist2

	for ind1=1:length(nsel)
		for ind2=1:length(nsel)
   			if ind1<ind2
				% A1=  zscore(subthresIndi(nsel(ind1),:));A2=  zscore(subthresIndi(nsel(ind2),:));
   				A1=  (subthresIndi(nsel(ind1),:));A2=  (subthresIndi(nsel(ind2),:));
   				[CC]=corrcoef(A1,A2,'Rows','complete');
 				% [c,lags]= xcorr(A1,A2,10,'Coeff');
 		  		allCmaxW(ind1,ind2)= CC(1,2);%max(abs(c));
 		   		allCmaxW(ind2,ind1)=allCmaxW(ind1,ind2);end
 		 		ROIdist2(ind1,ind2)= sqrt(sum((irm(:,nsel(ind1)).*sc -irm(:,nsel(ind2)).*sc).^2));
 		 		%ROIdist(ind1,ind2)=ROIdist(ind2,ind1);
			end
		end
   
	allCmaxW(allCmaxW>0.99)=NaN;
	allCmaxW(allCmaxW==0)=NaN;
   
	cwide= [cwide; allCmaxW(:)];
	rdist2= [rdist2; ROIdist2(:)];

	%% Calculate the lags of the spiking times
	subthresIndi= widefile.allresults.roaster;
	%subthresWide= widefile.allresults.roaster;

	clear allCmaxW  ROIdistW
	for ind1=1:length(nsel)
		for ind2=1:length(nsel)
	   		%if ind1<ind2
	  		A1=  zscore(subthresIndi(nsel(ind1),:));A2=  zscore(subthresIndi(nsel(ind2),:));
	 		
			if  length(find(A1>0)) >5 & length(find(A2>0)) >5
				%[c,lags]= xcorr(fastsmooth(A1,10,1,1),fastsmooth(A2,10,1,1),10,'Coeff');
				%allCmax(ind1,ind2)= max(abs(c));
	         		[c,lagsW]= xcorr(fastsmooth(A1,10,1,1),fastsmooth(A2,10,1,1),200,'Coeff');
	    			c1=c(191:211);
	    			[n1 n2]= max(abs(c1));
	     			
				if ind1~=ind2
	    				allCmaxW(ind1,ind2)= mean(c1);
	    				allCxW= [allCxW ; c];
	    				allCmW=[allCmW;c1(n2)];
	    				allRW= [allRW; sqrt(sum((irm(:,nsel(ind1))-irm(:,nsel(ind2))).^2))];
	     			else
	         			allCmaxW(ind1,ind2)= NaN;
	          			allCxW= [allCxW ; c.*NaN];
	            			allCmW=[allCmW;NaN];
	                		allRW= [allRW; NaN];
	   			end
	    
	 		else
				allCmaxW(ind1,ind2)=NaN;
			end
	  		
			% ROIdist(ind1,ind2)= sqrt(sum((irm(:,nsel(ind1))-irm(:,nsel(ind2))).^2));
		end
	end

	allCmaxW(allCmaxW>0.99)=NaN; %0.35 V1

	%   clear allCmaxW  
	% for ind1=1:length(nsel)
	% for ind2=1:length(nsel)
	%   A1=  zscore(subthresWide(nsel(ind1),:));A2=  zscore(subthresWide(nsel(ind2),:));
	%    if  length(find(A1>0)) >5 & length(find(A2>0)) >5
	%      [c,lags]= xcorr(fastsmooth(A1,10,1,1),fastsmooth(A2,10,1,1),10,'Coeff');
	%     allCmaxW(ind1,ind2)= max(abs(c)); else;  allCmaxW(ind1,ind2)=NaN;end
	% 
	% end;end
	%     allCmaxW(allCmaxW>0.6)=NaN;  
	%     
	    
	    
	cwideS= [cwideS; allCmaxW(:)];
	%cwideS= [cwideS; allCmaxW(:)];
  
end

rdist=rdist(~isnan(cindi));
cindi=cindi(~isnan(cindi));


rdist2=rdist2(~isnan(cwide));
cwide=cwide(~isnan(cwide));


% Plot the time lag for individual DMD
 figure('COlor','w','Position',[ 300 300 300 250],'Renderer', 'painters')
col=colormap(parula(7))
col=col(end:-1:1,:);
nj=0;
for id=100:100:600
    nj=nj+1;
plot(lags,nanmean(allCx(allR>id & allR<id+100,:)),'Color',col(nj,:));
hold on;
end; xlabel('Time lag'); ylabel('Corr')
legend( {num2str(150*0.32) ; num2str(250*0.32) ;num2str(350*0.32); num2str(450*0.32); num2str(550*0.32); num2str(650*0.32)})
title_string = 'Cross-correlation of spike-spike in Individual DMD';
title(title_string);

saveas(gcf, [save_fig_path 'Crosscorr lag\Jpeg Format\' title_string '.jpg']);
saveas(gcf, [save_fig_path 'Crosscorr lag\EPS Format\' title_string '.eps'], 'epsc');
saveas(gcf, [save_fig_path 'Crosscorr lag\SVG Format\' title_string '.svg']);


% Plot the time lag for widefield
 figure('COlor','w','Position',[ 300 300 300 250],'Renderer', 'painters')
col=colormap(parula(7))
col=col(end:-1:1,:);
nj=0;
for id=100:100:600
    nj=nj+1;
plot(lagsW,nanmean(allCxW(allRW>id & allRW<id+100,:)),'Color',col(nj,:));
hold on;
end; xlabel('Time lag'); ylabel('Corr')
legend( {num2str(150*0.32) ; num2str(250*0.32) ;num2str(350*0.32); num2str(450*0.32); num2str(550*0.32); num2str(650*0.32)})
title_string = 'Cross-correlation of spike-spike in Wide Field';
title(title_string);

saveas(gcf, [save_fig_path 'Crosscorr lag\Jpeg Format\' title_string '.jpg']);
saveas(gcf, [save_fig_path 'Crosscorr lag\EPS Format\' title_string '.eps'], 'epsc');
saveas(gcf, [save_fig_path 'Crosscorr lag\SVG Format\' title_string '.svg']);
