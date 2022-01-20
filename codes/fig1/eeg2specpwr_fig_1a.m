% Perform moving-window FFT on Raw EEG timeseries with Laplacian Referencing
% Adapted from Luis' 'NP_Step2_FastFourier_Periods_LuisCiria.m'
%
% update on 08/06/20: Bolun Chen
% 1. modified the code to automatically extract timeblock and latency for a
%    given physiological stage
% 2. applied moving-window FFT to generate spectrograms in every stages for all subjects 
% 3. applied Fourier Phase Randomization on the entire raw EEG data before applying moving-window FFT
%    results saved in '../SpecPow_PhaseRand/'
% 
% Updates on 08/07/20: Tested on pre-processed eeg data using Luis_Step1.m
%                      Block index overflows first occurs for
%                      '902_Test1_Task2';
%                      Also, noted that the pre-process data are not
%                      exactly the same as those copied from Luis folder
%                      (same trend & shape, smaller peaks in eeg
%                      timeseries)
%                     
% 
% 03/11/2021; Bolun Chen
% updated code computing spectral power from EEG signals
% 
% notes add on 05/03/2021: This matlab script contains all essential
% computation for checking consistency of fft, xcorr distribution and phase
% synchrony for Luis EEG data sets (Fig. 1 in the manuscript)
% Need a clean-up and better comments.
%

clear; close all; clc
%% add path to use eeglab
addpath '/Users/bolunchen/Documents/MATLAB/eeglab/'
eeglab
close
%% color maps
labelcolor=[0.9,    0.9,    0.9;        % R1
            0.68,   0.92,   1;          % WU
            1,      0.8,    0.8;        % EX
            0.68,   0.92,   1;          % CD
            1,      0.85,    0.5;       % T1        % 1,      0.9,    0.6;        % T1
            0.9,    0.9 ,   0.9;        % R2
            1,      0.85,    0.5];      % T2        % 1,      0.9,    0.6];       % T2

% colormap matrices of different freq-bands groups
histcolor_delta_grp = [0.8,     0,      0.2;        % delta-theta
                       0.7,     0.1,    0.1;        % delta-alpha
                       1,       0.2,    0;          % delta-sigma
                       1,       0.4,    0;          % delta-beta
                       1,       0.6,    0];         % delta-gamma
histcolor_theta_grp = [0,       1,      0;          % theta-alpha
                       0,       0.8,    0.2;        % theta-sigma
                       0,       0.6,    0.4;        % theta-beta
                       0,       0.4,    0.6];       % theta-gamma
histcolor_alpha_grp = [0,       0,      1;          % alpha-sigma
                       0,       0.3,    0.9;        % alpha-beta
                       0,       0.5,    0.8];       % alpha-gamma
histcolor_sigma_grp = [0.5,     0,      0.5;        % sigma-beta 
                       0.5,     0.3,    0.5];       % sigma-gamma
histcolor_beta_grp =  [1,       0.6,    0.7];       % beta-gamma

histcolor = [histcolor_delta_grp; histcolor_theta_grp; histcolor_alpha_grp; ...
             histcolor_sigma_grp; histcolor_beta_grp];        
        

%% parameters
FLAG_test = 1;                  % selecte Test 1 or Test 2 to analyze
FLAG_PhaseRand = 0;             % 0: raw eeg data;  1: perform phase randomization on eeg data

Fs = 1000;                      % Sampling frequency (Hz)
win = 2;                        % Window for FFT (sec)
overlap = 1;                    % Overlap time (sec)
win_len = win*Fs;               % # of data points in a window
mov_step = (win-overlap)*Fs;    % # of data in each moving step

num_freq = floor(win_len/2)+1;  % # of frequencies up to Nyquist freq (half of sampling rate)
                                % equals to (# data points in a window)/2 + 1 (zero freq)
min_freq = 0;                   % minimal frequency
max_freq = Fs/2;                % maximal frequency (Nyquist)
freq_vec = linspace(min_freq,max_freq,num_freq);    % physical frequencies (Hz)
dfreq = (max_freq-min_freq)/(num_freq-1);           % frequency spacing (Hz)

% six bands labels
freqbands = {'\delta', '\theta', '\alpha', '\sigma', '\beta', '\gamma'};
% frequency bands ranges (Hz)
delta_range = [0.5, 3.5]; 
theta_range = [4, 7.5];
alpha_range = [8, 11.5];
sigma_range = [12, 15.5];
beta_range = [16, 19.5];
gamma_range = [20, 24.5];    
% not distinguish gamma1 & gamma2: gamma1=[20,33], gamma2=[33.5,98.5], as in Rossella's work to account for high-freqcomponents of muscle activities;  

% freq-band boundary values
freq_range = [delta_range; theta_range; alpha_range; ...
              sigma_range; beta_range; gamma_range];          
% number of frequency bands
num_band = length(freq_range);  % with gamma, 6 bands; without gamma, 5 bands;

% freq-band boundary indices
freq_band_idx = floor(freq_range./dfreq) + 1;  


stages = {'Resting1', 'WarmUp', 'Exercise', ...
          'CoolDown', 'Task1', 'Resting2', 'Task2'};    % physiological stages
evnt_trigger = {'4','5','6','7','8','9','10'};          % event triggers for start of each stage
timeblk = [15, 10, 30, 10, 6, 15, 6].*60;               % duration of each stage in sec

chn = [1:9, 11:20, 22:31, 36];                          % selected channels out of 36 channels

%% Load database files

data_path = '../data/raw_eeg_data/';
data_filename_test1 = {dir([data_path,'*Test1*.set']).name};
data_filename_test2 = {dir([data_path,'*Test2*.set']).name};

if length(data_filename_test1) == length(data_filename_test2)
    num_subject = length(data_filename_test1);
else
    num_subject = min(length(data_filename_test1),length(data_filename_test2));
    disp('unequal data files for test 1 and test 2!!')
end
% NOTE: Only 18 subjects did both Test 1 & 2. 20 subjects did Test 2, among 
% them are sub-301 and sub-302 who didn't do Test1

subject_id = cell(1,num_subject);
for sub_idx = 1:num_subject
    marker = find(data_filename_test1{sub_idx}=='_');
    subject_id{sub_idx} = data_filename_test1{sub_idx}(1:marker(1)-1);
end

num_datafile = num_subject*2;   % num_sub*2, including 2 tests
num_stage = length(stages);     % # of stages
num_chn = length(chn);          % # of selected 30 channels

% create sub-directories to store spectrograms in each stage
switch FLAG_PhaseRand
    case 0
        parent_dir_name = './SpecPow/';
    case 1
        parent_dir_name = './SpecPow_PhaseRand/';
end
for stage_idx=1:num_stage
    save_dir = [parent_dir_name, stages{stage_idx}, '/'];
    mkdir(save_dir)
end

%% compute spectral power for six frequency bands
% Plot Fig.1a 

close all
% sub-figure position settings using stage durations
stage_dur = [15,10,30,10,6,15,6];  % duration of each stage in minutes
stage_gap = [7,0.5*ones(1,6)]; %ones(1,7);
% stage_gap = [5,0.5*ones(1,6)]; %ones(1,7);

sum_stage_dur = sum(stage_dur + stage_gap);    % total duration of 7 stages
sfig_width = stage_dur/sum_stage_dur;
sfig_hgap = stage_gap/sum_stage_dur;

left = cumsum(sfig_hgap) + cumsum([0, sfig_width(1:end-1)]); % - [0,0.03*ones(1,6)]
width = sfig_width;
% sub-figure position settings
% left = [0.08 0.21 0.30 0.55 0.64 0.698 0.828];
bottom = [0.85 0.71 0.57 0.43 0.29 0.15 0.1] - 0.075;
% width = [0.12 0.08 0.24 0.08 0.048 0.12 0.048];
height = ones(1,7)*0.125;
% height = ones(1,7)*0.13;

% stage labels
stage_labels = {'R1', 'WU', 'EX', 'CD', 'CT1', 'R2', 'CT2'}; 

% load spectral power data sets 
test = 'Test1';
root_dir = '../data/SpecPow/';
% output path
output_path = '../data/band_specpower_data/';


% for sub_idx = 1:num_subject
for sub_idx = 3;        % subject 1101
% for sub_idx = 4  % subject 1201
% for sub_idx = 10  % subject 1902
% for sub_idx = 12  % subject 201


    % subject ID
    subject = subject_id{sub_idx}
    
    % open new figure (with a fixed channel)
    figure(sub_idx)
    figtitle = strcat('sub-',subject,'-',test,'-C3');
    sgtitle(figtitle);  % subplot grid title

    for stage_idx = 1:num_stage
        sub_dir = strcat(stages{stage_idx},'/');

        data_path = [root_dir,sub_dir];
        data_filename = {dir([data_path,strcat(subject,'_',test,'*.mat')]).name} ;

        load([data_path,data_filename{1}])
        
        ch_idx = 8;  % channel C3
        chn_label = chn_labels{ch_idx};  % channel label
        % extract full spectral power time series at a given channel
        pwr_mat = squeeze(specpower(ch_idx,:,:))';      

        % time vector spectral power time series
        num_win = size(pwr_mat,2);      % # of window (time points for spec pwr)
        spec_time_vec = 1:num_win;      % time of moving windows (2sec)
        % define band spec power matrix
        band_pwr_mat = nan(num_band,num_win);   
        % define 14-point smoothed spec power matrix
        rel_band_pwr_smooth = nan(num_band,num_win);

        % compute band spec power
        for band_idx = 1:num_band
            % extract indices in the physical frequency vector for each freq-band
            freq_lims = freq_band_idx(band_idx,1):freq_band_idx(band_idx,2);

            % sum fourier modes' power over involved freqs in each band
%             band_pwr_mat(band_idx,:) = sum(pwr_mat(freq_lims,:),1);

            % avg fourier power over freqs in each band
            band_pwr_mat(band_idx,:) = mean(pwr_mat(freq_lims,:),1);
        end

        % compute relative band spec power, normalized to total sum of six bands
        rel_band_pwr_mat = band_pwr_mat./sum(band_pwr_mat,1);

        % smooth spectrogram with 14-point moving avg (for plot)
        smooth_span = 14;
        for band_idx = 1:num_band
            rel_band_pwr_smooth(band_idx,:) = ...
                smooth(rel_band_pwr_mat(band_idx,:),smooth_span);
        end

        for band_idx = 1:num_band
            pos = [left(stage_idx), bottom(band_idx), width(stage_idx), height(stage_idx)];
            subplot('Position',pos)   

            plot(rel_band_pwr_mat(band_idx,:))
%             plot(rel_band_pwr_mat(band_idx,:),'Color',[1 1 1]*0.6)

            hold on
            plot(rel_band_pwr_smooth(band_idx,:),'r','linew',1.2);
%             plot(rel_band_pwr_smooth(band_idx,:),'k','linew',1.2);


            xlim([0 num_win]); 
%             ylim([-0.05 1.05]); 
            [pwr_min,pwr_max] = bounds(rel_band_pwr_mat(band_idx,:));
            ylim([pwr_min-0.05 pwr_max+0.05])
            set(gca,'xtick',[],'ytick',[]);

            if(stage_idx==1)
%                 ylabel(strcat('$\tilde{S}(',freqbands{band_idx},')$'),'Interpreter','latex')
                ylabel(strcat('\boldmath$\widetilde{S}(',freqbands{band_idx},')$'),'Interpreter','latex')

%                 set(gca,'ytick',[0 1])
                set(gca,'ytick',[round(pwr_min,1) round(pwr_max,1)])

            end

            if(band_idx==1)
%                 xlabel(stage_labels{stage_idx})
                xlabel(strcat(num2str(stage_dur(stage_idx)),' min'))    
                set(gca,'XAxisLocation','top')
            end
            set(gca,'Color',labelcolor(stage_idx,:));   % set background color to distinguish different states;  to set them all grey, try: set(gca,'Color',[1,1,1]*0.94);                                   
        
        end
        
        % convert relative band power and smoothed matrices into tables
        rel_band_pwr_table = array2table(rel_band_pwr_mat','VariableNames',freqbands);
        rel_band_pwr_smooth_table = array2table(rel_band_pwr_smooth','VariableNames',freqbands);
        
        % save these tables as Excel spreadsheets
        exel_loc = '../graph_data_excel_sheets/';
        excel_name = 'fig_1a.xlsx';
        sheet_name = stages{stage_idx}; 
        
        writetable(rel_band_pwr_table, strcat(exel_loc,excel_name), ...
                   'Sheet',strcat(sheet_name,'_rel_band_pwr'))
        writetable(rel_band_pwr_smooth_table, strcat(exel_loc,excel_name), ...
                   'Sheet',strcat(sheet_name,'_rel_band_pwr_smooth'))

    end
      % save relative spec power plots at a fixed channel) as an .eps file
%       saveas(gcf,[output_path,figtitle,'_adapt_yaxis'],'epsc')
%       saveas(gcf,[output_path,figtitle],'epsc')
    saveas(gcf,'fig_1a','epsc')

end

%% fft function
function [freq_vec,data_fft,amp,pwr] = fft_func(time_vec,data_vec,sampling_rate,FLAG_plt)

    data_pnts = length(data_vec);           % # of data points
    data_fft = fft(data_vec)/data_pnts;     % scale fourier modes by # of data points
    
%     amp = 2*abs(data_fft);    % should exclude zero frequency component
    amp = abs(data_fft);
    amp(2:end) = 2*amp(2:end);              % double each frequency's component except for the zero frequency mode
    pwr = amp.^2;                           % power of each fourier mode 
    
    % Scale fourier modes' indices to physical freqencies in Hz
    % up to the Nyquist frequency (sampling_rate/2) 
    % # of fourier modes = # of data points
    min_freq = 0;
    max_freq = sampling_rate/2;
    num_freq = floor(data_pnts/2)+1;
    freq_vec = linspace(min_freq,max_freq,num_freq);
    
    % dfreq = (max_freq-min_freq)/(num_freq-1);
    % this is (sampling_rate/2)/floor(data_pnts/2), which is approximately, 
    % dfreq = sampling_rate/data_pnts = 1/signal_duration
 
    % visualize signals and fourier modes if FLAG_plt==1
    if FLAG_plt==1
       
        figure
        subplot(311)
        plot(time_vec,data_vec)
        xlabel('time (s)')
        ylabel('signal')

        subplot(312)
        data_ifft = ifft(data_fft)*data_pnts;
        plot(time_vec,data_ifft)
        xlabel('time (s)')
        ylabel('ifft[fourier modes]')

        subplot(313)
        plot(freq_vec,amp(1:num_freq),'-o',freq_vec,pwr(1:num_freq),'-s')
        xlim([min_freq min(30,max_freq)])
        xlabel('frequency (Hz)')
        ylabel('amplitude / power')
        legend({'amplitude','power'})
    end

end


%% define a function to trim EEG segements (removing rare spikes)
function [y,outlier_idx]= outlier_removal(x,m,method)
% input: data array x; multiple m (m*sigma to ignore)
% output: trimmed data array y; trimmed data indices idx

    y = x;
    [z,mean,~] = zscore(x);
    outlier_idx = find(abs(z)>m);
    
    switch method
        case 'mean'
            y(outlier_idx) = mean;
%             y(outlier_idx) = 0; 

        case 'interp'
            x_outs = [];
            if ~isempty(outlier_idx)        
                for out_idx=outlier_idx
                    if out_idx==1 
%                         x_out = interp1([out_idx+1,out_idx+2], x(out_idx+1:out_idx+2),out_idx,'linear');
                        x_out = interp1([out_idx+1,out_idx+2], x(out_idx+1:out_idx+2),out_idx,'spline');
                    elseif out_idx==length(x)
%                         x_out = interp1([out_idx-2,out_idx-1], x(out_idx-2:out_idx-1),out_idx,'linear');
                        x_out = interp1([out_idx-2,out_idx-1], x(out_idx-2:out_idx-1),out_idx,'spline');
                    else
%                         x_out = interp1([out_idx-1,out_idx+1], x(out_idx-1:2:out_idx+1),out_idx,'linear');
                        x_out = interp1([out_idx-1,out_idx+1], x(out_idx-1:2:out_idx+1),out_idx,'spline');
                    end
                    x_outs = [x_outs,x_out];
                end
                y(outlier_idx) = x_outs;
            end
    end
    
end
