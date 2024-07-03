function causal_features (subject_ids, path_input, path_output, working_dir)

%This function extracts the neural oscillations that are causally related to
%N-back performance using conditional independence-based causal sufficiency independent method

subject = sunject_ids;

for nid = 1:length(subject)
    
    disp('=============================================================')
    disp(num2str(nid))
    disp('=============================================================')
    pth_data=path_input; pth_subj='';
    %
    cd([pth_data,subject{nid},pth_subj]),pth=pwd;

    pth_data=path_input;
    
    cd([pth_data,subject{nid},'\3back']),pth=pwd;
    
    if isfile('power_baseline.mat')
        load(['perform.mat'])
        load(['power_baseline.mat'])
        load(['power_post.mat'])
        
        T2 = []; T2 = perform(:,2);
        
        [row, col] = find(T2<9);
        perform( row,:) = [];
        nb = 2;
        max_epoch=0;
        for cnd = 1:4
            if ~isempty (Rspns(nb+1).cond)
                if ~isempty(Rspns(nb+1).cond(cnd).v)
                    max_epoch (cnd) = max(Rspns(nb+1).cond(cnd).v(1,1,:));
                end
            end
        end
        max_2back = max(max_epoch);
        
        nb=3;
         
        ind_sorted = [];
        ind = [];
        
        for cnd = 1:4
            if ~isempty(Rspns(nb+1).cond(cnd).v)
                some_vector(:,:) = [squeeze(Rspns(nb+1).cond(cnd).v(1,1,:)- max_2back),...
                    repmat(cnd,(length(squeeze(Rspns(nb+1).cond(cnd).v(1,1,:)))),1)];
                ind = [ind; some_vector];
                clear some_vector
            end
        end
        ind_sorted = sortrows(ind,1);
        ind_sorted(row,:) =[];
        
        %
        % %  h = gscatter(ind_sorted(:,1),perform(:,2),ind_sorted(:,2),'rkgb','.',20,'on')
        % %
        % % legend(h(1:4),'TC','TNC', 'NTC','NTNC')
        %
        % perform =[ perform,ind_sorted(:,2)];
        %
        num = [1:size(perform,1)]';
        
        cd(working_dir);

        Electrodes = {'FP1';'FPZ';'FP2';'AF3';'AF4';'F7';'F5';'F3';'F1';'FZ';'F2';'F4';
            'F6';'F8';'FT7';'FC5';'FC3';'FC1';'FCZ';'FC2';'FC4';'FC6';'FT8';'T7';
            'C5';'C3';'C1';'CZ';'C2';'C4';'C6';'T8';'TP7';'CP5';'CP3';'CP1';'CPZ';
            'CP2';'CP4';'CP6';'TP8';'P7';'P5';'P3';'P1';'PZ';'P2';'P4';'P6';'P8';...
            'PO7';'PO5';'PO3';'POZ';'PO4';'PO6';'PO8';'O1';'OZ';'O2'};
        elec=[];
        for i=1:length(Electrodes) % find rejected electrodes
            ind=find(strcmp(Electrodes(i),  char(Rspns(5).elec)));
            elec = [elec ind];
        end
       
        index_all = size(perform,1);
        
        for epoch = 1:(index_all)
            y0(epoch,:) = (log10(reshape(Power_baseline {1}(:,:,(epoch)),1,[])));
            y2(epoch,:) = (log10(reshape(power_post {1}(:,:,(epoch)),1,[])));
        end
        for i = 1: size(y0,1)
            y0_norm = (y0-min(y0))./(max(y0)-min(y0));
            y2_norm = (y2-min(y2))./(max(y2)-min(y0));
        end

        count = length (Rspns(5).elec);
        num = [1:size(perform,1)]';

        [causesFoundHsic] =   mainSimulation(y0_norm, y2_norm ,perform(:,1));
        freq_sum = [];
        region_sum = {};
        for i =1: length(causesFoundHsic)
            freq_sum(i) = floor(((causesFoundHsic(i)-1)/count))+1;
            elect_ind = mod(causesFoundHsic(i),count);
            if (elect_ind ==0)
                elect_ind = count;
            end
            region_sum{i} = Rspns(5).elec{elect_ind};
        end
        
        cd([path_output, subject{nid}]),pth=pwd;

        save ('region_sum','region_sum')
        save ('freq_sum','freq_sum')
        
        clear y0_norm
        clear y2_norm
        clear y0
        clear y2
    end
    clear y0_norm
    clear y2_norm
    clear y0
    clear y2
end
end
