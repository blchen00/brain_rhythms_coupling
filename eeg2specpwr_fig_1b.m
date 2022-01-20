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

%% Load database files

data_path = '../data/raw_eeg_data/';

% extract subject id
subject = subject_id{sub_idx}
% select Test 1 or Test 2
test = strcat('Test',num2str(FLAG_test));
% determine .set data file's name
data_name = {dir([data_path,strcat(subject,'_',test,'*.set')]).name} ;
% load eeg data into a MATLAB structure 'EEG'
EEG = pop_loadset('filename',data_name,'filepath',data_path);

% exclude non-EEG channel locations
EEGchan_idx = find(~cellfun(@isempty,{EEG.chanlocs.X}))    
%%
cogtask_counter = 0;  % initialize a counter for concatenate cog-task 1 & 2 plots

%%  Plot Fig.1b 
% load spectral power data sets 
test = 'Test2';
root_dir = '../data/SpecPow/';
% output path
% output_path = '../band_specpower_data/';

stage_labels = {'R1', 'WU', 'EX', 'CD', 'CT1', 'R2', 'CT2'}; 

% for sub_idx = 1:num_subject
% for sub_idx = 3;      % sub-1101
% for sub_idx = 10;     % sub-1902
% for sub_idx = [4,12]  % sub-1201, sub-201
for sub_idx = 12;

    % subject ID
    subject = subject_id{sub_idx}
    
%     for stage_idx = 1:num_stage
    for stage_idx = 7;  %[1,3,[5,7]]; % rest, excercise, [cogtask1 cogtask2]
        stage_label = stages{stage_idx};
        sub_dir = strcat(stage_label,'/');

        data_path = [root_dir,sub_dir];
        data_filename = {dir([data_path,strcat(subject,'_',test,'*.mat')]).name} ;

        load([data_path,data_filename{1}])
        
        % time vector spectral power time series
        num_win = size(specpower,2);      % # of window (time points for spec pwr)
        spec_time_vec = 1:num_win;      % time of moving windows (2sec)
        
        % define matrices for spec powers
        sp_abs = nan(num_chn,num_band,num_win);
        sp_rel = nan(num_chn,num_band,num_win);
        sp_rel_smooth = nan(num_chn,num_band,num_win);
                
        for ch_idx = 1:num_chn
            chn_label = chn_labels{ch_idx};  % channel label

            % extract full spectral power time series at a given channel
            pwr_mat = squeeze(specpower(ch_idx,:,:))';      

            % define band spec power matrix
            band_pwr_mat = nan(num_band,num_win);   
            % define 14-point smoothed spec power matrix
            rel_band_pwr_smooth = nan(num_band,num_win);

            % compute band spec power
            for band_idx = 1:num_band
                % extract indices in the physical frequency vector for each freq-band
                freq_lims = freq_band_idx(band_idx,1):freq_band_idx(band_idx,2);
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
            
            % save spec power data of all channels during a stage
            sp_abs(ch_idx,:,:) = band_pwr_mat;
            sp_rel(ch_idx,:,:) = rel_band_pwr_mat;
            sp_rel_smooth(ch_idx,:,:) = rel_band_pwr_smooth;
        end
        
        % z-scoring sp_rel and sp_abs in a coarse-grained moving window      
        cg_win = 180;  % coarse grained window length in sec
        num_cg_win = floor(num_win/cg_win);     % number of coarse grained windows in a physiological stage

        cg_sp_abs = nan(num_chn,num_band,num_cg_win);
        cg_sp_rel = nan(num_chn,num_band,num_cg_win);

        for wi = 1:num_cg_win
            idx_range = 1+cg_win*(wi-1):cg_win*wi;
            cg_sp_abs(:,:,wi) = mean(sp_abs(:,:,idx_range),3);
            cg_sp_rel(:,:,wi) = mean(sp_rel(:,:,idx_range),3);
        end
        
        % calculate z-scores based on 3-min moving avg & std        
        mov_avg_sp_abs = movmean(sp_abs,cg_win,3);
        mov_std_sp_abs = movstd(sp_abs,cg_win,1,3);
        z_mov_sp_abs = (sp_abs-mov_avg_sp_abs)./mov_std_sp_abs;

        mov_avg_sp_rel = movmean(sp_rel,cg_win,3);
        mov_std_sp_rel = movstd(sp_rel,cg_win,1,3);
        z_mov_sp_rel = (sp_rel-mov_avg_sp_rel)./mov_std_sp_rel;
        
        smooth_span = 14;
        z_mov_sp_rel_smooth = movmean(z_mov_sp_rel,smooth_span,3);
        % smooth absolute spectral power
        sp_abs_smooth = nan(num_chn,num_band,num_win);
        for ch_idx = 1:num_chn
            for band_idx = 1:num_band
                sp_abs_smooth(ch_idx,band_idx,:) = ...
                    smooth(squeeze(sp_abs(ch_idx,band_idx,:)),smooth_span);
            end
        end
        
        % relative spec power plot
        ch_idx = 8;  % channel C3
        
        figure
        tiledlayout(num_band,1,'TileSpacing','tight','Padding','compact')
        for band_idx = 1:num_band
            nexttile
            % plot relative spec power & its smooths
            plot(spec_time_vec,squeeze(sp_rel(ch_idx,band_idx,:)))
%             plot(spec_time_vec,squeeze(sp_rel(ch_idx,band_idx,:)),'color',[1 1 1]*0.6)
            hold on
            plot(spec_time_vec,squeeze(sp_rel_smooth(ch_idx,band_idx,:)),'r','linew',1.2)

            % plot moving avg, z-score & smoothed z-score of relative spec power in 3min window
            %{
            plot(spec_time_vec,squeeze(mov_avg_sp_rel(ch_idx,band_idx,:)))
            plot(spec_time_vec,squeeze(z_mov_sp_rel(ch_idx,band_idx,:)),'color',[1 1 1]*0.6)
            hold on
            plot(spec_time_vec,squeeze(z_mov_sp_rel_smooth(ch_idx,band_idx,:)),'linew',2)
            legend({'SP-rel','movmean','z-score'})
            %}       
            xlim([0 num_win]); xticklabels([])            
            % xlim([(num_win-cg_win)/2 (num_win+cg_win)/2]); % ylim([-0.05 1.05])
            [pwr_min,pwr_max] = bounds(squeeze(sp_rel_smooth(ch_idx,band_idx,:)));
%             [pwr_min,pwr_max] = bounds(squeeze(sp_rel(ch_idx,band_idx,:)));

            ylim([pwr_min-0.1 pwr_max+0.1])
%             ylim([pwr_min-0.05 pwr_max+0.05])

%             set(gca,'ytick',[round(pwr_min,1) round(pwr_max,1)])
            set(gca,'ytick',[])
            if(band_idx==num_band)
                xlabel('time (s)')
                xticklabels('auto')
            end
          
%             ylabel(strcat('$\tilde{S}_{',freqbands{band_idx},'}$'),'Interpreter','latex')
%             ylabel(strcat('$\tilde{S}(',freqbands{band_idx},')$'),'Interpreter','latex')
%             ylabel(strcat('\boldmath$\widetilde{S}_{',freqbands{band_idx},'}$'),'Interpreter','latex')

            set(gca,'Color',labelcolor(stage_idx,:))
        end        
        set(gcf,'Position',[616,845,318,172]);
%         saveas(gcf,strcat('sub_',subject,'_',test,'_',stage_label,'_SPrel'),'epsc')   
        
        % plot selected epochs in Exercise 
        
        if (stage_idx ==3)
           time2plt = 400:1600;
           
           figure
            tiledlayout(num_band,1,'TileSpacing','compact','Padding','compact')
            for band_idx = 1:num_band
                nexttile
                % plot relative spec power & its smooths
                plot(time2plt,squeeze(sp_rel(ch_idx,band_idx,time2plt)))
                hold on
                plot(time2plt,squeeze(sp_rel_smooth(ch_idx,band_idx,time2plt)),'r','linew',1.2)     
                xlim([time2plt(1) time2plt(end)]);
                xticks(linspace(time2plt(1),time2plt(end),6))
                xticklabels([])            
                % xlim([(num_win-cg_win)/2 (num_win+cg_win)/2]); % ylim([-0.05 1.05])
                [pwr_min,pwr_max] = bounds(squeeze(sp_rel_smooth(ch_idx,band_idx,time2plt)));
                ylim([pwr_min-0.05 pwr_max+0.05])
                set(gca,'ytick',[round(pwr_min,1) round(pwr_max,1)])
%                 set(gca,'ytick',[])
                if(band_idx==num_band)
                    xlabel('time (s)')
                    xticklabels('auto')
                end

                ylabel(strcat('$\tilde{S}_{',freqbands{band_idx},'}$'),'Interpreter','latex')
            %             ylabel(strcat('$\tilde{S}(',freqbands{band_idx},')$'),'Interpreter','latex')

                set(gca,'Color',labelcolor(stage_idx,:))
            end        
            %         set(gcf,'Position',[617,226,560,791]);
            set(gcf,'Position',[616,350,317,420]);
%             saveas(gcf,strcat('sub_',subject,'_',test,'_',stage_label,'_SPrel'),'epsc')       
        end
        
        
        % concatenate cogtask 1 & 2
        if (stage_idx==5)
            cogtask_counter = 1;
            cogtask_sp_rel = []; cogtask_sp_rel_smooth = []; spec_tvec =[];
            
            sp_rel_tmp = squeeze(sp_rel(ch_idx,:,:));
            sp_rel_smooth_tmp = squeeze(sp_rel_smooth(ch_idx,:,:));
            tvec_tmp = spec_time_vec;
            
            cogtask_sp_rel = sp_rel_tmp;
            cogtask_sp_rel_smooth = sp_rel_smooth_tmp;
            spec_tvec = tvec_tmp;
        end
            
        if (stage_idx==7 && cogtask_counter==1)
            sp_rel_tmp = squeeze(sp_rel(ch_idx,:,:));
            sp_rel_smooth_tmp = squeeze(sp_rel_smooth(ch_idx,:,:));
            tvec_tmp = spec_time_vec;
            
            cogtask_sp_rel = [cogtask_sp_rel, sp_rel_tmp(:,2:end)];
            cogtask_sp_rel_smooth = [cogtask_sp_rel_smooth, sp_rel_smooth_tmp(:,2:end)];
            
            spec_tvec = [spec_tvec, tvec_tmp(2:end)+spec_tvec(end)];
            
            % plot concatenated cog-task1 & 2 data
            figure
            tiledlayout(num_band,1,'TileSpacing','tight','Padding','compact')
            for band_idx = 1:num_band
                nexttile
                % plot relative spec power & its smooths
                plot(spec_tvec(1:end-1),cogtask_sp_rel(band_idx,1:end-1))
                hold on
                plot(spec_tvec(1:end-1),cogtask_sp_rel_smooth(band_idx,1:end-1),'r','linew',1.2)
                
                xlim([0 spec_tvec(end)]);
                xticks(linspace(0,spec_tvec(end),10))
                xticklabels([])            
                % xlim([(num_win-cg_win)/2 (num_win+cg_win)/2]); % ylim([-0.05 1.05])
                [pwr_min,pwr_max] = bounds(cogtask_sp_rel_smooth(band_idx,:));
                ylim([pwr_min-0.05 pwr_max+0.05])
    %             set(gca,'ytick',[round(pwr_min,1) round(pwr_max,1)])
                set(gca,'ytick',[])
                if(band_idx==num_band)
                    xlabel('time (s)')
                    xticklabels('auto')
                end

%                 ylabel(strcat('$\tilde{S}_{',freqbands{band_idx},'}$'),'Interpreter','latex')
                ylabel(strcat('\boldmath$\widetilde{S}_{',freqbands{band_idx},'}$'),'Interpreter','latex')

                set(gca,'Color',labelcolor(stage_idx,:))
                xline(0.5*spec_tvec(end),'k:')
            end        

            set(gcf,'Position',[616,845,318,172]);

%             saveas(gcf,strcat('sub_',subject,'_',test,'_','Task1-2','_SPrel'),'epsc')   
                 
            cogtask_sp_rel = []; cogtask_sp_rel_smooth = []; spec_tvec =[];
        end


        % convert relative band power and smoothed matrices into tables
        rel_band_pwr_table = array2table(squeeze(sp_rel(ch_idx,:,:))','VariableNames',freqbands);
        rel_band_pwr_smooth_table = array2table(squeeze(sp_rel_smooth(ch_idx,:,:))','VariableNames',freqbands);
        
        % save these tables as Excel spreadsheets
        exel_loc = '../graph_data_excel_sheets/';
        excel_name = 'fig_1b.xlsx';
        sheet_name = stages{stage_idx}; 
        file_name = strcat(exel_loc,excel_name);
        
        writetable(rel_band_pwr_table, file_name, ...
                   'Sheet',strcat(sheet_name,'_rel_band_pwr'))
        writetable(rel_band_pwr_smooth_table, file_name, ...
                   'Sheet',strcat(sheet_name,'_rel_band_pwr_smooth'))
    end

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
