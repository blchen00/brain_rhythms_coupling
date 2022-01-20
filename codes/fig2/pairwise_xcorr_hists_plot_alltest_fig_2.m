% Histogram Profiles & Smooth avg for freqband pairs in all stages
% Input:  pairwise correaltion matrices for all pairs
%         Need names for those matrice .mat files (around line 80)
% Output: One matlab figure (.fig & .eps) 
%         It should be saved as full screen (see the last section)
%
% NOTE!! Change the saved figure's name according to loaded .mat file!! 
% 07/18/20: Bolun Chen
% 
% update on 08/04/20: plotted correlation histogram profiles with 1000 shuffles
%                     used correlation data in '../../Stat_Tests_SP_shuffles/
%                     figure saved in '../../Figures/shuffles/'
% update on 08/19/20: Saved a copy for generating Fig.2 in
%                     /Final_figures_Bolun/codes/
%
clear; clc; close all
%% Basic Settings & Parameters
% set input & output paths
input_path = '../data/correlations/C3/';    % where input .mat locates
% output_path = '../raw_figs/';                 % where to save figures

% names of physiological stages
stages = {'Resting1', 'WarmUp', 'Exercise',  ...
          'CoolDown', 'Task1', 'Resting2', 'Task2'}; 
% names of frequency bands
freqbands = {'\delta', '\theta', '\alpha', '\sigma', '\beta', '\gamma'};

num_stage = length(stages);             % # of stages
num_band = length(freqbands);           % # of freq-band
num_pair = num_band*(num_band-1)/2;     % # of freq-band pairs

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
            1,      0.9,    0.6;        % T1
            0.9,    0.9 ,   0.9;        % R2
            1,      0.9,    0.6];       % T2

% histogram settings
bw = 0.05;                  % bin width
bmin = -1;  bmax = 1;       % bin limits
bin_edges = bmin:bw:bmax;   % bin edegs
smooth_span = 5; %4;            % step-size for smooth bincounts

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
         
%% Plot Correlation Histogram Profiles in Tiledlayout
fig_handle = figure('Position', get(0, 'Screensize'));
%{
tileplots = tiledlayout(num_stage+1,num_pair+1);
tileplots.TileSpacing = 'none';
tileplots.Padding = 'none';   
%}

panelplots = panel();
panelplots.pack('h', {5/num_pair, 7/num_pair, 3/num_pair})
panelplots(1).pack(num_stage,5)
panelplots(2).pack(num_stage,7)
panelplots(3).pack(num_stage,3)

%%
for s_idx=1:num_stage
    % pairwise correlation matrix .mat filename
    pairwise_xcorr_mat = strcat('pairwise_xcorr_chn_C3_sub_all_test_all_stage_', ...
                                stages{s_idx}, '.mat');

    % load pairwise correlation matrix .mat file from the designated directory
    load([input_path, pairwise_xcorr_mat]);
    

    % create table for each physiological state, starting with bin_edges
    hist_table = array2table(bin_edges', 'VariableNames',{'bin_edges'});
    hist_smooth_table = array2table(bin_edges', 'VariableNames',{'bin_edges'});
    % loop through every pair of freq-bands
    for p_idx = 1:num_pair
        % correlations between a pair of freq-bands
        C = corr_pairwise(:,p_idx);    
        % histogram(C, bin_edges, 'Normalization', 'pdf')

        bin_counts = histc(C,bin_edges);                        % should use histcounts instead!
        bin_counts = bin_counts./max(bin_counts);               % rescale by max bin_counts
        smooth_bin_counts = smooth(bin_counts,smooth_span);     % smooth the rescaled bin_counts
        
        % save bin_counts, smooth_bin_counts to excel sheets
        hist_table = addvars(hist_table, bin_counts, 'NewVariableNames',pairlabel{p_idx});
        hist_smooth_table = addvars(hist_smooth_table, smooth_bin_counts, 'NewVariableNames',pairlabel{p_idx});

        % save these tables as Excel spreadsheets
        exel_loc = '../graph_data_excel_sheets/';
        excel_name = 'fig_2.xlsx';
        sheet_name = stages{s_idx}; 
        file_name = strcat(exel_loc,excel_name);
        
        writetable(hist_table, file_name, 'Sheet',strcat(sheet_name,'_histogram'))
        writetable(hist_smooth_table, file_name, 'Sheet',strcat(sheet_name,'_histogram_smooth'))

        % locat tiledlayout coordinates
        %{
        row = s_idx;  col = p_idx;
        tile_idx = col + row*(num_pair+1) + 1;
        nexttile(tile_idx)
        %}
        
        % locat panel plots coordinates
        row = s_idx; 
        if p_idx<=5
            group = 1;
            col = p_idx;
        elseif p_idx<=12
            group = 2;
            col = p_idx-5;
        else 
            group = 3;
            col = p_idx-12;
        end
        panelplots(group,row,col).select();
        
        % histogram profile 
        hold on
        bar(bin_edges,bin_counts, 0.5, ...
            'EdgeColor',histcolor(p_idx,:), 'FaceColor',histcolor(p_idx,:));
        % smooth avg
        plot(bin_edges,smooth_bin_counts,'LineWidth',3,'Color','k');
        
        % add vertical dashed line at zero correlation
        xline(0,'k--','LineWidth',2)
%         xline(0.5,'k:','LineWidth',2.5)
%         xline(-0.5,'k:','LineWidth',2.5)
        
        % axis settings
        axis([-1.05 1.05 0 1]); 
        axis square
        xticks([]); yticks([]);
        box on; ax = gca; ax.LineWidth = 1;
        
        % add xticks & yticks to subplots
        %{
        set(gca,'XTickLabel',[],'YTickLabel',[])
        if s_idx==num_stage
            xticks([-1 0 1]); 
            ax = gca; ax.FontSize = 22;
        end
        if p_idx==1
            yticks([0 1]); 
            ax = gca; ax.FontSize = 22;
        end
        %}
        
    end

end

panelplots.de.margin = 6;
panelplots(1).marginleft = 0;
panelplots(2).marginleft = 18;
panelplots(3).marginleft = 18;

%% Adding labels in plots
%{
% plot frequency-band pairs labels
for p_idx = 1:num_pair
    tile_idx = p_idx + 1;  
    nexttile(tile_idx)
    
    axis([-1 1 0 1]); axis square; axis off

    text (-0.6,0,pairlabel{p_idx},'LineWidth',1, 'FontSize',25,...
        'EdgeColor',histcolor(p_idx,:), 'BackgroundColor','White');
end

% plot physiological states labels
for s_idx = 1:num_stage
    tile_idx = s_idx*(num_pair+1) + 1;
    nexttile(tile_idx)
    
    axis([-1 1 0 1]); axis square; axis off
    
    text(0,0.5,statelabel{s_idx},'LineWidth',1, 'FontSize',20,...
        'EdgeColor','Black', 'BackgroundColor',labelcolor(s_idx,:), ...
        'Rotation',90);
end
%}

%% Save figures
% should be saved as full screen, otherwise small plots are scrambled
fig_name = 'hist_profiles_chn_C3_sub_all_test_all_stage_all';
% savefig(fig_handle, [output_path, fig_name]);
% saveas(fig_handle,[output_path, fig_name],'epsc')