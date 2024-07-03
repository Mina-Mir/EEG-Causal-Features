function power_calculation(subject_ids, path_input, path_output,working_dir)

% This function calculates the pre- and post-stimuls powers from 4 frequncy
% bands.


TIME_cut=0;                %= set to 1 to accept epoch with reponse time less than 3sec

random_cond  =0;           % random_cond=0;  to concatenate epoches with TIME_cut condition sum([1:rt1, 1:rt2])==5sec

Sgmnt_rnd_lng=5000;        % length of random segmantation window
signal_length=5000;        % specify  signal length (ms)

srate    = 1000;
Electrod = {'FP1';'FPZ';'FP2';'AF3';'AF4';'F7';'F5';'F3';'F1';'FZ';'F2';'F4';
    'F6';'F8';'FT7';'FC5';'FC3';'FC1';'FCZ';'FC2';'FC4';'FC6';'FT8';'T7';
    'C5';'C3';'C1';'CZ';'C2';'C4';'C6';'T8';'TP7';'CP5';'CP3';'CP1';'CPZ';
    'CP2';'CP4';'CP6';'TP8';'P7';'P5';'P3';'P1';'PZ';'P2';'P4';'P6';'P8';...
    'PO7';'PO5';'PO3';'POZ';'PO4';'PO6';'PO8';'CB1';'O1';'OZ';'O2';'CB2'};

subject = subject_ids;

for nid=1:length(subject)
    
    pth_data=path_input; pth_subj='';
    
    cd([pth_data,subject{nid},pth_subj]),pth=pwd;
  
    Fs = 1000;
    disp('=============================================================')
    disp(num2str(nid))
    disp('=============================================================')
    EEG=[];
    
    cd(working_dir') %Folder address of the saved script
    
    cd(path_output),pth=pwd;
    load Rspns_c
    load perform
    cd(working_dir) %Folder address of the saved script
    
    T2 = []; T2 = perform(:,2);
    R2=[]; R2  = ref2mastoid(Rspns_c,Rspns(5).elec);
    [row, col] = find(T2<9);
    T2( row,:) = [];
    R2(:,:,row) = [];
    
    chan60 = [(1:57),(59:61)];
    
    chan=Rspns(5).elec;
    
    %  Excluding CB1, CB2
    
    ind = find(strcmp(Electrod(62), chan));
    if  ~isempty(ind), chan(ind) = [];end
    ind = find(strcmp(Electrod(58), chan));
    if  ~isempty(ind),  chan(ind) = [];end
    
    
    time_win = [1:4500];
    R22=[]; R22 = R2(:,time_win,:);
    epoch_size = size(R22,3);
    
    chan_size = size(R22,1);
    
    if (epoch_size ~=1)
        
        R_erp = [];
        R_erp(:,:,1)  = squeeze(mean (R22,3));
        R = []; R = zeros(chan_size,length(time_win),epoch_size); % subtracting the ERP
        for epoch = 1:epoch_size
            
            R(:,:,epoch) = R22(:,:,epoch) - R_erp;
        end
    else
        R=[]; R = R22;
    end
    
    % wavelet parameters
    min_freq = 1;
    max_freq = 50;
    num_frex = 100;
    srate = 1000;
    time_pnts = [-1400:701];
    time_baseline = [-700 -100];
    [~,baselineidx(1)]=min(abs(time_pnts-time_baseline(1)));
    [~,baselineidx(2)]=min(abs(time_pnts-time_baseline(2)));
    
    time_base = [baselineidx(1):baselineidx(2)];
    % other wavelet parameters
    
    frequencies = linspace(min_freq,max_freq,num_frex);
    time = -1:1/srate:1;
    half_of_wavelet_size = (length(time)-1)/2;
    
    % FFT parameters (use next-power-of-2)
    n_wavelet     = length(time);
    n_data        = length(time_pnts);
    n_convolution = n_wavelet+n_data-1;
    n_conv_pow2   = pow2(nextpow2(n_convolution));
    wavelet_cycles= 5;
    
    % FFT of data
    freq_range_ind = {7:13; 14:23; 24:58; 60:100};
    
    RT_m = 700;%round(T2(1,ind));
    
    TF_sum = zeros(chan_size,4,1400+RT_m);
    for epoch = 1:epoch_size
        tf_data = zeros(chan_size,length(frequencies),length(time_pnts));
        for chan=1:chan_size
            
            fft_data = fft(squeeze(R(chan,:,epoch)),n_conv_pow2);
            
            % initialize output time-frequency data
            
            for fi=1:length(frequencies)
                
                % create wavelet and get its FFT
                wavelet = (pi*frequencies(fi)*sqrt(pi))^-.5 * exp(2*1i*pi*frequencies(fi).*time) .* exp(-time.^2./(2*( wavelet_cycles /(2*pi*frequencies(fi)))^2))/frequencies(fi);
                fft_wavelet = fft(wavelet,n_conv_pow2);
                
                % run convolution
                convolution_result_fft = ifft(fft_wavelet.*fft_data,n_conv_pow2);
                convolution_result_fft = convolution_result_fft(1:n_convolution); % note: here we remove the extra points from the power-of-2 FFT
                convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
                
                % put power data into time-frequency matrix
                tf_data(chan,fi,:) = (abs(convolution_result_fft).^2)*fi;
            end
        end
        for i =1:4
            tf_data_avf_freq(:,i,:) =squeeze( mean ( tf_data(:,freq_range_ind{i},:),2));
            
            
        end
        
        tf_data_avf_freq_rt_post(:,:,:) = tf_data_avf_freq(:,:,1501:1400+RT_m);
        
        
        for freq = 1:4
            for chan = 1: chan_size
                Power_baseline{1}(chan,freq,epoch) = squeeze(mean (tf_data_avf_freq(chan,freq,time_base),3));
                power_post{1} (chan,freq,epoch) =  squeeze(mean (tf_data_avf_freq_rt_post(chan,freq,:),3));
                
            end
        end
    end
    
    clear tf_data_avf_freq
    clear tf_data_avf_freq_rt_post
    
    cd(path_output),pth=pwd;
    save ('Power_baseline','Power_baseline')
    save ('power_post','power_post')
    
end
end