% plot fig 1c 

clear; close all; clc
%% parameters
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
pos_scatter_color = [0 0 0];
neg_scatter_color = [1 1 0.8];

scatter_size = 150;
xcorr_win_length = 30; 
%% load pairwise correlation data file for a given subject
data_path = '../data/correlations/pairwise_xcorr/C3/subjects/test1/';
subject_id = '102';

filenames = {dir(data_path).name};
mat_names = filenames(contains(filenames,subject_id,'IgnoreCase',true));

% load data & make plots for Resting1 & Exercise
for s_idx = [1,3]
    selected_stage = stages{s_idx}
    % load .mat file for the given stage
    mat_name = mat_names{contains(mat_names, selected_stage, 'IgnoreCase', true)};
    load([data_path,mat_name])
    
    plt_background_color = labelcolor(s_idx,:);
    
    figure
    fig = tiledlayout(3,1);
    for p_idx = [5,6,15]
%         pairlabel{p_idx}
        % extract correlation for a given pair
        xcorr = corr_pairwise(:,p_idx);
        pos_idx = find(xcorr>0);
        xcorr_pos = xcorr(pos_idx);
        neg_idx = find(xcorr<=0);
        xcorr_neg = xcorr(neg_idx);
    
        nexttile
        plot(xcorr,'k-','LineWidth',1.5)
        hold on
        scatter(pos_idx,xcorr_pos,scatter_size,...
                'MarkerEdgeColor',[0 0 0],...
                'MarkerFaceColor',pos_scatter_color,'LineWidth',1.5)
        scatter(neg_idx,xcorr_neg,scatter_size,...
                'MarkerEdgeColor',[0 0 0],...
                'MarkerFaceColor',neg_scatter_color, 'LineWidth',1.5)
        yline(0,'--','LineWidth',1)
        ylim([-1.2 1.2])
        yticks([-1 0 1])
        yticks([])
        bands = split(pairlabel{p_idx},'-');
        ylabel(strcat('$C(',bands{1},'$-$',bands{2},')$'),'Interpreter','latex')
        xticks(linspace(0,length(xcorr),4))
        xticklabels(linspace(0,length(xcorr),4)*xcorr_win_length)
        set(gca,'FontSize',20,'Color',plt_background_color)
    end
    fig.TileSpacing = 'compact';
    fig.Padding = 'compact';
    sgtitle(strcat('subject-',subject_id," test1 chnC3 xcorr ",selected_stage))

    % convert correlation array to a table
    %{
    time_vec = (1:size(corr_pairwise,1))'*xcorr_win_length;
    xcorr_table = array2table([time_vec, corr_pairwise(:,[5,6,15])], ...
                               "VariableNames",{'time(sec)',pairlabel{[5,6,15]}});

    % save these tables as Excel spreadsheets
    exel_loc = '../graph_data_excel_sheets/';
    excel_name = 'fig_1c.xlsx';
    sheet_name = strcat(selected_stage,'_pairwise_correlation'); 
    file_name = strcat(exel_loc,excel_name);
    
    writetable(xcorr_table, file_name, 'Sheet',sheet_name)
    %}
end

%% load data & make plots for Task1 & Task2 combined
corr_pairwise_combined = [];
for s_idx = [5,7]
    selected_stage = stages{s_idx}
    % load .mat file for the given stage
    mat_name = mat_names{contains(mat_names, selected_stage, 'IgnoreCase', true)};
    load([data_path,mat_name])
    
    corr_pairwise_combined = [corr_pairwise_combined; corr_pairwise];
end

plt_background_color = labelcolor(s_idx,:);
    
figure
fig = tiledlayout(3,1);
for p_idx = [5,6,15]
%     pairlabel{p_idx}

    xcorr = corr_pairwise_combined(:,p_idx);
    pos_idx = find(xcorr>0);
    xcorr_pos = xcorr(pos_idx);
    neg_idx = find(xcorr<=0);
    xcorr_neg = xcorr(neg_idx);

    nexttile
    plot(xcorr,'k-','LineWidth',1.5)
    hold on
    scatter(pos_idx,xcorr_pos,scatter_size,...
            'MarkerEdgeColor',[0 0 0],...
            'MarkerFaceColor',pos_scatter_color,'LineWidth',1.5)
    scatter(neg_idx,xcorr_neg,scatter_size,...
            'MarkerEdgeColor',[0 0 0],...
            'MarkerFaceColor',neg_scatter_color, 'LineWidth',1.5)
    yline(0,'--','LineWidth',1)
    xline(0.5*length(xcorr),':','LineWidth',1)
    ylim([-1.2 1.2])
    yticks([-1 0 1])
    bands = split(pairlabel{p_idx},'-');
    ylabel(strcat('$C(',bands{1},'$-$',bands{2},')$'),'Interpreter','latex')
    xticks(linspace(0,length(xcorr),4))
    xticklabels(linspace(0,length(xcorr),4)*xcorr_win_length)
    set(gca,'FontSize',20,'Color',plt_background_color)
end
fig.TileSpacing = 'compact';
fig.Padding = 'compact';

sgtitle(strcat('subject-',subject_id," test1 chnC3 xcorr Tasks1-2"))

% convert correlation array to a table
%{
time_vec = (1:size(corr_pairwise_combined,1))'*xcorr_win_length;
xcorr_table = array2table([time_vec, corr_pairwise_combined(:,[5,6,15])], ...
                           "VariableNames",{'time(sec)',pairlabel{[5,6,15]}});

% save these tables as Excel spreadsheets
exel_loc = '../graph_data_excel_sheets/';
excel_name = 'fig_1c.xlsx';
sheet_name = 'Task1-2_pairwise_correlation'; 
file_name = strcat(exel_loc,excel_name);

writetable(xcorr_table, file_name, 'Sheet',sheet_name)
%}

%% generate raw figures with the right size for creating the figure panel (fig. 1c)
close all
scatter_size = 100;

% load data & make plots for Resting1 & Exercise
for s_idx = [1,3]
    selected_stage = stages{s_idx}
    % load .mat file for the given stage
    mat_name = mat_names{contains(mat_names, selected_stage, 'IgnoreCase', true)};
    load([data_path,mat_name])
    
    plt_background_color = labelcolor(s_idx,:);

    for p_idx = [5,6,15]
        figure

        % extract correlation for a given pair
        xcorr = corr_pairwise(:,p_idx);
        pos_idx = find(xcorr>0);
        xcorr_pos = xcorr(pos_idx);
        neg_idx = find(xcorr<=0);
        xcorr_neg = xcorr(neg_idx);
        
        % plot correlation values for a given pair
        plot(xcorr,'k-','LineWidth',1.5)
        hold on
        scatter(pos_idx,xcorr_pos,scatter_size,...
                'MarkerEdgeColor',[0 0 0],...
                'MarkerFaceColor',pos_scatter_color,'LineWidth',1)
        scatter(neg_idx,xcorr_neg,scatter_size,...
                'MarkerEdgeColor',[0 0 0],...
                'MarkerFaceColor',neg_scatter_color, 'LineWidth',1)
        yline(0,'--','LineWidth',1)
        ylim([-1.2 1.2])
        xlim([0 length(xcorr)+1])
        set(gca,'Color',plt_background_color,'XTick',[],'YTick',[])
        set(gcf,'Units','Inches','Position',[0, 0, 5, 5*0.4])
    end

end


% load data & make plots for  Task1 & Task2 combined
corr_pairwise_combined = [];
for s_idx = [5,7]
    selected_stage = stages{s_idx}
    % load .mat file for the given stage
    mat_name = mat_names{contains(mat_names, selected_stage, 'IgnoreCase', true)};
    load([data_path,mat_name])
    
    corr_pairwise_combined = [corr_pairwise_combined; corr_pairwise];
end

plt_background_color = labelcolor(s_idx,:);
    
for p_idx = [5,6,15]
    figure

    xcorr = corr_pairwise_combined(:,p_idx);
    pos_idx = find(xcorr>0);
    xcorr_pos = xcorr(pos_idx);
    neg_idx = find(xcorr<=0);
    xcorr_neg = xcorr(neg_idx);

    plot(xcorr,'k-','LineWidth',1.5)
    hold on
    scatter(pos_idx,xcorr_pos,scatter_size,...
            'MarkerEdgeColor',[0 0 0],...
            'MarkerFaceColor',pos_scatter_color,'LineWidth',1)
    scatter(neg_idx,xcorr_neg,scatter_size,...
            'MarkerEdgeColor',[0 0 0],...
            'MarkerFaceColor',neg_scatter_color, 'LineWidth',1)
    yline(0,'--','LineWidth',1)
    xline(0.5*length(xcorr),':','LineWidth',1)
    ylim([-1.2 1.2])
    xlim([0 length(xcorr)+1])
    set(gca,'Color',plt_background_color,'XTick',[],'YTick',[])
    set(gcf,'Units','Inches','Position',[0, 0, 5, 5*0.4])
end



%% calculate histograms (popualtion average) of the above pairwise correlations
%{
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

panelplots = panel();
panelplots.pack('h', {5/num_pair, 7/num_pair, 3/num_pair})
panelplots(1).pack(num_stage,5)
panelplots(2).pack(num_stage,7)
panelplots(3).pack(num_stage,3)

%% generate pannel plots
% data location
input_path = '../data/correlations/pairwise_xcorr/C3/';

for s_idx=1:num_stage
    % pairwise correlation matrix .mat filename
    pairwise_xcorr_mat = strcat('pairwise_xcorr_chn_C3_sub_all_test_1_stage_', ...
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
        excel_name = 'fig_1c.xlsx';
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
%}
