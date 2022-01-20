% Generate bar charts of positive & negative correaltion ratios 
% Inputs:   input_path, subjecct_input_path, output_path
%           subject_IDs & chn_labels information
%           pairwise correlation histogram data
% 
% Outputs:  corr_ratios;    % cell array for corr-ratio of each subject
%           corr_ratios_grp;% 3-dim matrix for grp avg & std of corr-ratio
%           bar-chart figure; 
% 
% 07/21/20: Bolun Chen
%
% updates on 08/19/20: plotted C3 (test 1 & 2) barplots 
%
clear; clc; close all
%% Basic Settings & Parameters
% set input & output paths
input_path = '../data/xcorr_spnorm/';       % where input .mat locates
subject_input_path = '../data/pairwise_xcorr/C3/subjects/';    % pariwise xcorr including test 1
output_path = '../raw_figs/group_avg/';     % where to save figures

% names of physiological stages
stages = {'Resting1', 'WarmUp', 'Exercise',  ...
          'CoolDown', 'Task1', 'Resting2', 'Task2'}; 
% names of frequency bands
freqbands = {'\delta', '\theta', '\alpha', '\sigma', '\beta', '\gamma'};

num_stage = length(stages);             % # of stages
num_band = length(freqbands);           % # of freq-band
num_pair = num_band*(num_band-1)/2;     % # of freq-band pairs

% load subject IDs & channel labels
xcorr_data_prefix = 'xcorr_spnorm_alltest_allsub_allchn_gamma_';

xcorr_data = strcat(xcorr_data_prefix, stages{1},'.mat');
load([input_path, xcorr_data], 'subject_IDs','chn_labels');
num_sub = length(subject_IDs);
num_channel = length(chn_labels);

% labels for states & freqband pairs
statelabel = {'R1';'WU';'EX';'CD';'CT1';'R2';'CT2'};
pairlabel = cell(1,num_pair);
count = 0;
for i=1:num_band-1
    for j=i+1:num_band
        pair = strcat(freqbands{i},  '-', freqbands{j});
        count = count+1;
        pairlabel{count} = pair;
    end
end       

% colormap matrix of state-labels
labelcolor=[0.9,    0.9,    0.9;        % R1
            0.68,   0.92,   1;          % WU
            1,      0.8,    0.8;        % EX
            0.68,   0.92,   1;          % CD
            1,      0.85,    0.5;       % T1        % 1,      0.9,    0.6;        % T1
            0.9,    0.9 ,   0.9;        % R2
            1,      0.85,    0.5];      % T2        % 1,      0.9,    0.6];       % T2

% histogram settings
bw = 0.05;                  % bin width
bmin = -1;  bmax = 1;       % bin limits
bin_edges = bmin:bw:bmax;   % bin edegs
smooth_span = 5; %4;            % step-size for smooth bincounts

% critical values for postive & negative correlations
C_pos_val = 0.5;    C_neg_val = -0.5;
% indices of ranges for postive & negative correlation
bin_edges_pos = find(bin_edges>=C_pos_val);
bin_edges_neg = find(bin_edges<=C_neg_val);


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
         
%% Plot Group Average & Std in Barcharts

fig_handle = figure('Position', [100 100 400 1000]);

panelplots = panel();
panelplots.pack('v', num_stage)

corr_ratios = cell(num_stage, num_pair);        % all subjects corr_ratio vectors
corr_ratios_grp = nan(4,num_pair, num_stage);   % group avg & std for pos & neg ratios

for s_idx=1:num_stage
    
    for p_idx = 1:num_pair
        
        corr_ratios_sub = nan(num_sub,2);
        % plot each subject's avg correlation profile on top of group avg
        for sub_idx = 1:num_sub

            % load correaltion data .mat for one subject 
            subject_xcorr_mat = strcat('pairwise_xcorr_chn_C3_sub_', subject_IDs{sub_idx}, ...
                                        '_test_all_stage_', stages{s_idx}, '.mat');

            load([subject_input_path, subject_xcorr_mat]);
            
            C_sub = corr_pairwise(:,p_idx);                                     % note the group .mat corr_pairwise has been replaced
            bin_counts_sub = histc(C_sub,bin_edges);                            % should use histcounts instead!
            bin_counts_sub = bin_counts_sub./max(bin_counts_sub);               % rescale by max bin_counts
            smooth_bin_counts_sub = smooth(bin_counts_sub,smooth_span);         % smooth the rescaled bin_counts
            
            total_bin_area_sub = sum(smooth_bin_counts_sub)*bw;                 % total area under smoothed bin_counts
            pos_bin_area_sub = sum(smooth_bin_counts_sub(bin_edges_pos))*bw;    % area of postive bin_counts
            neg_bin_area_sub = sum(smooth_bin_counts_sub(bin_edges_neg))*bw;    % area of negative bin_counts
            pos_corr_ratio_sub = pos_bin_area_sub./total_bin_area_sub;           % ratio of positive correlations over entire histogram
            neg_corr_ratio_sub = neg_bin_area_sub./total_bin_area_sub;           % ratio of negative correlations over entirw histogram
                     
            corr_ratios_sub(sub_idx,:) = [(-1).*neg_corr_ratio_sub, pos_corr_ratio_sub];
        end
        
        % save each subject's corr profile ratio vectors to a cell array
        corr_ratios{s_idx, p_idx} = corr_ratios_sub;
        
        % group avg & std of corr profile ratios
        corr_ratios_mean = mean(corr_ratios_sub);   % [neg_corr_ratio, pos_corr_ratio]
        corr_ratios_std = std(corr_ratios_sub);
        % save grp avg & std to a cell array
        corr_ratios_grp(:, p_idx, s_idx) = [corr_ratios_mean'; corr_ratios_std']; 
    end
    
    % plot bar charts for each stage & freq-band pair    
    panelplots(s_idx).select();
    
    hold on; box on;
    neg_corr_ratio_avg = corr_ratios_grp(1,:,s_idx);
    pos_corr_ratio_avg = corr_ratios_grp(2,:,s_idx);
    
    bar(neg_corr_ratio_avg, 'FaceColor', [1 1 0.8])
    bar(pos_corr_ratio_avg, 'FaceColor', labelcolor(s_idx,:))
    
    % error bars
    neg_corr_ratio_std = corr_ratios_grp(3,:,s_idx);
    pos_corr_ratio_std = corr_ratios_grp(4,:,s_idx);

    er = errorbar(neg_corr_ratio_avg, neg_corr_ratio_std); 
    er.Color = 'black';    
    er.LineStyle = 'none';  

    er = errorbar(pos_corr_ratio_avg, pos_corr_ratio_std);    
    er.Color = 'black';                            
    er.LineStyle = 'none';  
    
    % axis settings
    axis([0.5 15.5 -1 1]);
    xticks(1:num_pair); yticks(-1:0.5:1);
    ax = gca; ax.LineWidth = 1;
    ax.XTickLabel = []; ax.FontSize = 10;

    % add freq-band pair labels inside bar charts
    pos_bar_tips = pos_corr_ratio_avg+pos_corr_ratio_std; 
    neg_bar_tips = neg_corr_ratio_avg-neg_corr_ratio_std; 

    for p_idx = 1:num_pair
        if p_idx<=5
            text(p_idx-0.45,pos_bar_tips(p_idx) + 0.1, ...
                 pairlabel{p_idx},'FontSize',15)
        else
            text(p_idx-0.5,neg_bar_tips(p_idx) - 0.075 - 0.135*mod(p_idx-1,2), ... 
                 pairlabel{p_idx},'FontSize',15)
        end
    end

    % add auxiliary lines as guides to the eyes
    yline(0.5,'k:');  yline(-0.5,'k:');
    xline(5.5,'k','LineWidth',0.7);
%     xline(9.5,'k','LineWidth',0.7);
    xline(12.5,'k','LineWidth',0.7);
%     xline(14.5,'k','LineWidth',0.7);

    % add physiological stage labels
    %{
    text (1,0.5, statelabel{s_idx}, 'LineWidth',0.5, 'FontSize',8,...
    'EdgeColor','Black','BackgroundColor',labelcolor(s_idx,:));
    %}
    
    % remove xticks
    set(gca,'XTick',[])
end

%% adjusting margins between pannels
panelplots.de.margin = 5;
panelplots.fontsize = 16;

%% Save figures
fig_name = 'barchart_C3_allsub_alltest_group_avg';
savefig(fig_handle, [output_path, fig_name]);
saveas(fig_handle,[output_path, fig_name],'epsc')

%% Save each subject's corr profile ratio and group avg&std
%{
% corr_ratios = cell(num_stage, num_pair);        % all subjects corr_ratio vectors
% corr_ratios_grp = nan(4,num_pair, num_stage);   % group avg & std for pos & neg ratios

barchart_cellarray_mat = strcat('barchart_cellarray_C3_allsub_alltest');
save_dir = '../data/';                        
save([save_dir, barchart_cellarray_mat], 'corr_ratios','corr_ratios_grp');
%}


%% plot barchart for each subject
%{
% for sub_idx = 1:num_sub

for sub_idx = 9   % subject '902'; Luis used subject '102' (sub_idx=1)

    figure('Position', [100 100 400 1000]);
    
    panelplots = panel();
    panelplots.pack('v', num_stage)

    for s_idx = 1:num_stage
%         subplot(num_stage,1,s_idx)
        
        panelplots(s_idx).select();
        hold on; box on;

        neg_corr_ratio = [];
        pos_corr_ratio = [];
        for p_idx = 1:num_pair
            neg_corr_ratio = [neg_corr_ratio, corr_ratios{s_idx,p_idx}(sub_idx,1)];
            pos_corr_ratio = [pos_corr_ratio, corr_ratios{s_idx,p_idx}(sub_idx,2)];
        end

        bar(neg_corr_ratio, 'FaceColor', [1 1 0.8])
        bar(pos_corr_ratio, 'FaceColor', labelcolor(s_idx,:))

         % axis settings
        axis([0.5 15.5 -1 1]);
        xticks(1:num_pair); yticks(-1:0.5:1);
        ax = gca; ax.LineWidth = 1;
        ax.XTickLabel = []; ax.FontSize = 10;
%         axis tight
        set(gca,'XTick',[])
        
        % add freq-band pair labels inside bar charts
        pos_bar_tips = pos_corr_ratio; 
        neg_bar_tips = neg_corr_ratio; 

        for p_idx = 1:num_pair
            if p_idx<=5
                text(p_idx-0.45,pos_bar_tips(p_idx) + 0.15, ...
                     pairlabel{p_idx},'FontSize',15)
            else
                text(p_idx-0.5,neg_bar_tips(p_idx)- 0.125 - 0.15*mod(p_idx-1,2), ...
                     pairlabel{p_idx},'FontSize',15)
            end
        end
        
        % add auxiliary lines as guides to the eyes
        yline(0.5,'k:');  yline(-0.5,'k:');
        xline(5.5,'k','LineWidth',0.7);
%         xline(9.5,'k','LineWidth',0.7);
        xline(12.5,'k','LineWidth',0.7);
%         xline(14.5,'k','LineWidth',0.7);

    end
%% adjusting margins    
    panelplots.de.margin = 5;
    panelplots.fontsize = 16;

 %% save figure   
    saveas(gcf,[output_path, 'barchart_C3_sub_', subject_IDs{sub_idx}],'epsc')
end
%}