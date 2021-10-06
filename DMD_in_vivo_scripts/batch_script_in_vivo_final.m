clear all

%% CHANGE paths to data and where the processed files needs to be saved
savepath1=['\\engnas.bu.edu\research\eng_research_handata\EricLowet\DMD\comp_wide_indi\'];

path1='\\engnas.bu.edu\research\eng_research_handata\Yangyang Wang\Rotation-DMD invivo data processing\July 4 processing by folder\organized by individual & large mask\';
cd(path1)

%%
main_folder = {'A'; 'B';'C';'D'}

time_win=300; % in 2ms bins (=600ms), window at the start of the trial and end of the trial  to estimate phottobleaching
camera_noise=158.8; %  camera noise variance for a given pixel

for id=1:length(main_folder)  %%%% MAIN LOOP
    cd([path1 main_folder{id}])
    
    subF=dir;
    for id2=1:length(subF)
        if length(subF(id2).name )>3 & length(strfind(subF(id2).name,'txt'))==0
            cd([path1 main_folder{id} '\' subF(id2).name   ])
            
            if length(dir('traces*'))>0   %%%%%%%%%%%%%
                
                ro= dir('ROI*');
                load(ro(1).name);
                clear allresults
                sesname= dir('traces*');
                allbleach=[];  allbleach2=[];  allbleach3=[];  allbleach4=[];  allbleach5=[];  allbleach6=[];
                for  idd=1:length(sesname)
                    load(sesname(idd).name)
                    clear neuronS
                    for xd=1:length(ROIs)
                        neuronS(xd)= (camera_noise)./length(find(ROIs{xd}));
                    end
                    %
                    result=spike_detect_SNR_final(traces,neuronS);
                    
                    if idd==1
                        allresults=result;
                    else
                        f = fieldnames( allresults);
                        for i = 1:length(f)
                            allresults.(f{i}) =[   allresults.(f{i}),result.(f{i})];
                        end
                        allresults.trial_l   =[     allresults.trial_l , result.trial_l];
                    end
                    [ble] = photobleach_estimation_final(traces, time_win);
                    ble(:,1)=(1-ble(:,1));
                    allbleach=[allbleach; ble(:,1)];
                    
                end
                allresults.roi=ROIs;
                allresults.bleach=allbleach;
          
                if ~isempty(  strfind(subF(id2).name, 'wide'))
                    allresults.type= 'wide';
                    save([savepath1  subF(id2).name '_allresults_wide' '.mat'], 'allresults')
                else
                    allresults.type= 'indi';
                    save([savepath1  subF(id2).name '_allresults_indi' '.mat'], 'allresults')
                end
                % save([savepath1  subF(id2).name '_allresults.mat'], 'allresults')
                
                
            else %%%%%%%%%%%%%%%%%%%%%%%%
                subF2=dir;
                for id3=1:length(subF2)   % 
                    if length(subF2(id3).name )>3 & length(strfind(subF2(id3).name,'txt'))==0
                        cd([path1 main_folder{id} '\'  subF(id2).name '\' subF2(id3).name])
                        if length(dir('traces*'))>0
                            ro= dir('ROI*');
                            load(ro(1).name);
                            clear allresults
                            sesname= dir('traces*');
                            allbleach=[];  allbleach2=[];  allbleach3=[];  allbleach4=[];  allbleach5=[];  allbleach6=[];
                            for  idd=1:length(sesname)
                                load(sesname(idd).name)
                                
                                clear neuronS
                                for xd=1:length(ROIs)
                                    neuronS(xd)= (camera_noise)./length(find(ROIs{xd}));
                                end
                                %
                                result=spike_detect_SNR_final(traces,neuronS);
                                if idd==1
                                    allresults=result;
                                else
                                    f = fieldnames( allresults);
                                    for i = 1:length(f)
                                        allresults.(f{i}) =[   allresults.(f{i}),result.(f{i})];
                                    end
                                    allresults.trial_l   =[     allresults.trial_l , result.trial_l];
                                end
                                [ble] = photobleach_estimation_final(traces, time_win);
                                ble(:,1)=(1-ble(:,1));
                                
                                allbleach=[allbleach; ble(:,1)];
                       
                                
                            end
                            allresults.bleach=allbleach;
                            allresults.roi=ROIs;
                          
                            if ~isempty(  strfind(subF2(id3).name, 'wide'))
                                allresults.type= 'wide';
                                save([ savepath1  subF(id2).name '_' subF2(id3).name '_allresults_wide' '.mat'], 'allresults');
                            else
                                allresults.type= 'indi';
                                save([ savepath1  subF(id2).name '_' subF2(id3).name '_allresults_indi' '.mat'], 'allresults');
                            end
                            
                        else
                            
                            subF3=dir;
                            for id4=1:length(subF3)
                                if length(subF3(id4).name )>3 & length(strfind(subF3(id4).name,'txt'))==0
                                    cd([path1 main_folder{id} '\'  subF(id2).name '\' subF2(id3).name '\' subF3(id4).name])
                                    if length(dir('traces*'))>0
                                        ro= dir('ROI*');
                                        load(ro(1).name);
                                        clear allresults
                                        sesname= dir('traces*');
                                        allbleach=[];  allbleach2=[];  allbleach3=[];  allbleach4=[];  allbleach5=[];  allbleach6=[];
                                        for  idd=1:length(sesname)
                                            load(sesname(idd).name)
                                            clear neuronS
                                            for xd=1:length(ROIs)
                                                neuronS(xd)= (camera_noise)./length(find(ROIs{xd}));
                                            end
                                            %
                                            result=spike_detect_SNR_final(traces,neuronS);
                                            if idd==1
                                                allresults=result;
                                            else
                                                f = fieldnames( allresults);
                                                for i = 1:length(f)
                                                    allresults.(f{i}) =[   allresults.(f{i}),result.(f{i})];
                                                end
                                                allresults.trial_l   =[     allresults.trial_l , result.trial_l];
                                            end
                                            [ble] = photobleach_estimation_final(traces, time_win);
                                            ble(:,1)=(1-ble(:,1));
                                            allbleach=[allbleach; ble(:,1)];
                                         
                                            
                                        end
                                        allresults.bleach=allbleach;
                                          allresults.roi=ROIs;
                                        
                                        
                                        if ~isempty(  strfind(subF3(id4).name, 'wide'))
                                            allresults.type= 'wide';
                                            save([ savepath1  subF(id2).name '_' subF2(id3).name '_' subF3(id4).name '_allresults_wide' '.mat'], 'allresults');
                                        else
                                            allresults.type= 'indi';
                                            save([ savepath1  subF(id2).name '_' subF2(id3).name '_' subF3(id4).name '_allresults_indi' '.mat'], 'allresults');
                                        end
                                        
                                    end  % 3rd level
                                end
                            end
                            
                        end  % 2nd level
                    end
                end
            end %1st level
        end
    end
    
    
end %main loop
