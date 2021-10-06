
% In vivo directory
%cd('\\engnas.bu.edu\research\eng_research_handata\EricLowet\DMD\invivoDMD\')
%cd('\\engnas.bu.edu\research\eng_research_handata\EricLowet\DMD\comp_wide_indi\')

% In vitro directory
cd('\\ad\eng\research\eng_research_handata\Pierre Fabris\DMD Project\All In Vitro Analysis\');

ses=dir('*.mat');

%%% search for indi
clear findwide
for id=1:length(ses)
    if   strfind(ses(id).name, 'Wide')>0
        findwide(id)=0;
    else
        findwide(id)=1;
    end
    
end
indiloc=find(findwide);
%%% search for fitting indi

%% seelct matching ROI
%indiB=[];wideB=[]; indiSNR=[]; wideSNR=[];
cindi=[];cwide=[]; rdist=[];allR=[];
cindiS=[];cwideS=[]; allCx=[];allCm=[]; allPairs=[];
for id=1:length(indiloc)
    id
  F=  ses(indiloc(id)).name;
%if any([~isempty(strfind(ses(indiloc(id)).name, '612777'))  ~isempty(strfind(ses(indiloc(id)).name, '613200')) ]==1 ) % V1 selectin
 %  if any([~isempty(strfind(F, '2020.05.29'))   ~isempty(strfind(F, '608442'))   ~isempty(strfind(F, '602086'))    ~isempty(strfind(F, '602088'))  ]==1 ) % hippo, many cells
  %  if any([~isempty(strfind(F, '608445'))    ]==1 ) % striate
   if any([~isempty(strfind(F, 'Culture 8_Wide Field'))   ]==1 )

    %   widefile=load(ses(wideloc(id)).name);
    indifile=load(ses(indiloc(id)).name);

    clear irm
    for id2= 1:length(indifile.allresults.roi)
        [ x y]=find(indifile.allresults.roi{id2});
        irm(:,id2)= round(mean([x , y]));end
    ROIs=indifile.allresults.roi;

    
%    cindi= [cindi; allCmax(:)];
%cwide= [cwide; allCmaxW(:)];
%rdist= [rdist; ROIdist(:)];

    %%
    subthresIndi= indifile.allresults.roaster;

    
if 0  % try out .... ignore that part
       subthresIndi= indifile.allresults.trace_ws;
 
 subthresIndi(isnan(subthresIndi))=0;
 
 
lfp=[];[chas, timed, tr]=size( subthresIndi)
for trial=1:tr
    trial
    lfp.trial{trial}= (double( subthresIndi(:,1:end,trial)));
     lfp.trial{trial}= bsxfun(@minus,lfp.trial{trial}, nanmean(lfp.trial{trial},2));
    lfp.time{trial} = (1:size( lfp.trial{trial},2))./500;
end
for ch=1:chas
   lfp.label{ch}= [ 'A' num2str(ch)];
end

 cfg1=[];
        cfg1.channel= 'all';
        cfg1.bpfilter='yes';
            cfg1.bpfreq=[  10 12];%  
     [lfpA] = ft_preprocessing(cfg1,lfp);
     for id=1:size(lfp.trial,2)
      for id2=1:size(lfp.trial{1},1)
       lfp.trial{id}(id2,:)  =   lfp.trial{id}(id2,:) - lfpA.trial{id}(id2,:) ;
      end;end
   cfg1=[];
        cfg1.channel= 'all';
        cfg1.bpfilter='yes';
            cfg1.bpfreq=[  3 9];%  
            cfg.hilbert='angle'
     [lfpA] = ft_preprocessing(cfg1,lfp);
     for id=1:size(lfp.trial,2)
      for id2=1:size(lfp.trial{1},1)
       lfp.trial{id}(id2,:)  =  lfp.trial{id}(id2,:);% real(lfpA.trial{id}(id2,:)) ;
      end;end
  
   subthresIndi=lfp.trial{1};
end
 %  subthresWide= widefile.allresults.roaster;
 nsel=1:size(subthresIndi,1);
    clear allCmax  
for ind1=1:length(nsel)
    ind1
for ind2=1:length(nsel)
   %   if ind1<ind2
  A1=  zscore(subthresIndi(nsel(ind1),:));A2=  zscore(subthresIndi(nsel(ind2),:));
 if  length(find(A1>0)) >1 & length(find(A2>0)) >1  %minimum spike number used 20

     
corwin=400; % time lag
         [c,lags]= xcorr(fastsmooth(A1,3,1,1),fastsmooth(A2,3,1,1),corwin,'Coeff');
    c1=c(corwin-4:corwin+5);
    [n1 n2]= max(abs(c1));
     if ind1~=ind2
    allCmax(ind1,ind2)= mean(c1);
    allCx= [allCx ; c];
    allCm=[allCm;c1(n2)];
    allPairs=[allPairs; [ nsel(ind1) nsel(ind2)]];
    allR= [allR; sqrt(sum((irm(:,nsel(ind1))-irm(:,nsel(ind2))).^2))];
     else
         allCmax(ind1,ind2)= NaN;
          allCx= [allCx ; c.*NaN];
            allCm=[allCm;NaN];
                allR= [allR; NaN];  allPairs=[allPairs; [nsel(ind1) nsel(ind2)]];
   end
    
 else;  allCmax(ind1,ind2)=NaN;end
  % ROIdist(ind1,ind2)= sqrt(sum((irm(:,nsel(ind1))-irm(:,nsel(ind2))).^2));
end;end
    allCmax(allCmax>0.99)=NaN;
cindiS= [cindiS; allCmax(:)];
%cwideS= [cwideS; allCmaxW(:)];

 end
  
end
%rdist2=rdist;
rdist=rdist(~isnan(cindi));
cindi=cindi(~isnan(cindi));



%rdist2=rdist2(~isnan(cindiS))
cindiS=cindiS(~isnan(cindiS))

% PLOT


%% Average image
load('\\engnas.bu.edu\research\eng_research_handata\Pierre Fabris\DMD Project\v2 script Trace Extraction\In vitro\Culture 8\Wide Field\average_somArchon_1.mat');

figure('COlor','w')
imagesc(averageFrame)
colormap(gray);hold on,
for ind=1:size(irm,2)
%plot(irm(2,ind).*1,irm(1,ind).*1,'.r'); hold on;

text(irm(2,ind),irm(1,ind), sprintf('%d', ind), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle','COlor', [ 1 1 1]);
end
set(gca,'Xticklabel', [],'Yticklabel',[])    






%%%%%%%%%%%%%%%%%%


CN=averageFrame;
CN=CN-min(CN(:));
CN=CN./max(CN(:));
clear CN3;
CN3(:,:,1)= CN;
CN3(:,:,2)= CN;
CN3(:,:,3)= CN;
As= zeros(size(CN3,1), size(CN3,2),3);
As(:,:,1)=1;As(:,:,2)=0.;
ref=11; % Reference ROI
tag=[ 1 4 13] ;% target ROI
figure('COlor','w'),subplot(1,1,1),imshow(CN3);hold on,
for ref2=1:max(allPairs(:))
t=find(allPairs(:,1)==ref  & allPairs(:,2)==ref2);
if ~isempty(t)
sel=[ref  ref2];
v1=abs(allCm(t));
As= zeros(size(CN3,1), size(CN3,2),3);
As(:,:,1)=v1*5;As(:,:,2)=v1*5;
h=imshow(As);%hold off;
AF= ROIs{sel(1)};AF2= ROIs{sel(2)};set(h,'AlphaData', (AF2)./1.5)
end
end
plot(irm(2,ref).*1,irm(1,ref).*1,'ow','Markersize',23); 
plot(irm(2,ref).*1,irm(1,ref).*1,'ow','Markersize',29); 
%plot(irm(2,ref).*1,irm(1,ref).*1,'ow','Markersize',33);
for ind=[ ref tag]
text(irm(2,ind),irm(1,ind), sprintf('%d', ind), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle','COlor', [ 1 1 1]);
end
for j=1:length(tag)
plot(irm(2,tag(j)).*1,irm(1,tag(j)).*1,'ow','Markersize',23); 
end
set(gca,'Xticklabel', [],'Yticklabel',[])
figure('COlor','w'), 
for j=1:length(tag)
subplot(length(tag),1,j)
plot(lags*2,allCx(find(allPairs(:,1)==ref & allPairs(:,2)==tag(j)),:),'k');axis tight;ylim([ -0.02 0.11])
title([ num2str(ref) '    '   num2str(tag(j)) ]);xlabel('lag ms');ylabel('Corr')
end


% 
% figure,
% for ind=1:77
% plot(nanmean(allCx(find(allPairs(:,1)==ind ),:),1)+ind.*0.02)
% hold on,
% end
% 
% figure,nt=0;
% for ref2=1:77;
% ref=62;
% t=find(allPairs(:,1)==ref  & allPairs(:,2)==ref2)
% if ~isempty(t)
%   vx=  allCm(t);
%   if vx>0.04;nt=nt+.04;
%     plot(lags,allCx(find(allPairs(:,1)==ref & allPairs(:,2)==ref2),:)+nt)
%     hold on;
%   end
% end
% end;%xlim([-10 10])
%     
