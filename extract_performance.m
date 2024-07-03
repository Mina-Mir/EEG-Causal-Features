function extract_performance (subject_ids, path_input, path_output, working_dir)

% this function calculates the N-back performance using the marks (keypad response and event types) on EEG
% data after ICA decomposition


subject = subject_ids;

TIME_cut=0;                %= set to 1 to accept epoch with reponse time less than 3sec
nid=1;                          %= set to 0 to accept epcoh with reponse time more than 3sec


random_cond  =0;           % random_cond=0;  to concatenate epoches with TIME_cut condition sum([1:rt1, 1:rt2])==5sec
% random_cond=1;  to create a signal with length =signal_length
% random_cond==2  to create a signal with length =signal_length,but non -overlapping segmantation

Sgmnt_rnd_lng=5000;        % length of random segmantation window
signal_length=5000;        % specify  signal length (ms)


srate    = 1000;
Electrod = {'FP1';'FPZ';'FP2';'AF3';'AF4';'F7';'F5';'F3';'F1';'FZ';'F2';'F4';
    'F6';'F8';'FT7';'FC5';'FC3';'FC1';'FCZ';'FC2';'FC4';'FC6';'FT8';'T7';
    'C5';'C3';'C1';'CZ';'C2';'C4';'C6';'T8';'TP7';'CP5';'CP3';'CP1';'CPZ';
    'CP2';'CP4';'CP6';'TP8';'P7';'P5';'P3';'P1';'PZ';'P2';'P4';'P6';'P8';...
    'PO7';'PO5';'PO3';'POZ';'PO4';'PO6';'PO8';'CB1';'O1';'OZ';'O2';'CB2'};
%=====================================================================================================================
%=====================================================================================================================
%%
for nid=1:length(subject)
    
    clc
    disp('=============================================================')
    disp(num2str(nid))
    disp('=============================================================')
    
    pth_data=path_input;

    cd (pth_data)


    cd([pth_data,subject{nid}]),pth=pwd;

    dd=dir('*_RC.set');

    if ~isempty(dd)

        EEG = pop_loadset('filename',dd.name); 
        
        L=ones(4,5); nev=length(EEG.event); EEG.event(nev+1).type='END';EEG.event(nev+2).type='END';
        Rspns={};
        
        %%
        cd(working_dir)
        % =========== create Response matrix ======================================
        %  accpeted reponse for RT<TIME and not-switched reponses
        
        if TIME_cut==1; TIME=3000; else TIME=4500; end
        t = 4500;
        L=ones(4,5); nev=length(EEG.event); EEG.event(nev+1).type='END';EEG.event(nev+2).type='END';
        Rspns={};
        Rspns_c = [];
        perform = [];
        count = 1;
        nb=3;

        for ev=1:length(EEG.event)
            
            disp([num2str(nid),'   ',num2str(nb),'    ',num2str(ev)])
            if strcmp([num2str(nb),'*'],EEG.event(ev).type(1:2))
                if strcmp('10',{EEG.event(ev+1).type})
                    if strcmp('keypad1',{EEG.event(ev+2).type})
                        RT=EEG.event(ev+2).latency-EEG.event(ev+1).latency;
                        if (RT<TIME)&&(~strcmp('keypad2',{EEG.event(ev+3).type}))
                            perform (count,1) = 1;
                            perform (count,2) = RT;
                            Rspns_c (:,1:t,count) = EEG.data(:,:,EEG.event(ev).epoch);
                            count = count+1;
                      
                        end
                    elseif strcmp('keypad2',{EEG.event(ev+2).type})
                        RT=EEG.event(ev+2).latency-EEG.event(ev+1).latency;
                        if (RT<TIME)&&(~strcmp('keypad1',{EEG.event(ev+3).type}))
                            perform (count,1) = 0;
                            perform (count,2) = RT;
                            Rspns_c (:,1:t,count) = EEG.data(:,:,EEG.event(ev).epoch);
                            count = count+1;
                        end
                    end
                end
                if strcmp('20',{EEG.event(ev+1).type})
                    if strcmp('keypad1',{EEG.event(ev+2).type})
                        RT=EEG.event(ev+2).latency-EEG.event(ev+1).latency;
                        if (RT<TIME)&&(~strcmp('keypad2',{EEG.event(ev+3).type}))
                            perform (count,1) = 0;
                            perform (count,2) = RT;
                            Rspns_c (:,1:t,count) = EEG.data(:,:,EEG.event(ev).epoch);
                            count = count+1;
                            
                        end
                    elseif strcmp('keypad2',{EEG.event(ev+2).type})
                        RT=EEG.event(ev+2).latency-EEG.event(ev+1).latency;
                        if (RT<TIME)&&(~strcmp('keypad1',{EEG.event(ev+3).type}))
                            perform (count,1) = 1;
                            perform (count,2) = RT;
                            Rspns_c (:,1:t,count) = EEG.data(:,:,EEG.event(ev).epoch);
                            count = count+1;
                            
                        end
                    end
                end
   
                
            end

        end  % for ev
        
        
    end
    mkdir([path_output,subject{nid}])
    cd([path_output, subject{nid}]),pth=pwd;
    save ('perform', 'perform');
    save('Rspns_c','Rspns_c')
end




