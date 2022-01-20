% Generate network graphs from barcharts avg correlation profiles 
% Inputs:   input_path,  output_path
%           barcharts data: corr_ratios_grp, corr_ratios;
% 
% Outputs:  networks graphs 
% 
% 12/31/20: Bolun Chen
% update: 1/12/21   1) change color-coding to be consistent with previous
%                   figures
%
clear; clc; close all
%% Basic Settings & Parameters
% set input & output paths
input_path = '../data/barcharts_data/';       % where input .mat locates
output_path = '../raw_figs/networks/';     % where to save figures

% names of physiological stages
stages = {'Resting1', 'WarmUp', 'Exercise',  ...
          'CoolDown', 'Task1', 'Resting2', 'Task2'}; 
% channels' labels
chn_labels = {'Fp1', 'Fp2', 'C3', 'C4', 'O1', 'O2'};
      
% names of frequency bands
freqbands = {'\delta', '\theta', '\alpha', '\sigma', '\beta', '\gamma'};

num_stage = length(stages);             % # of stages
num_band = length(freqbands);           % # of freq-band
num_pair = num_band*(num_band-1)/2;     % # of freq-band pairs
num_channel = length(chn_labels);       % # of channels

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
         
% color maps for markers
markercolor = [0        1       1; 
               0        1       0; 
               1        0.5     0.5;
               1        1       0;  
               0.9      0.5 	0;
               1        0       0; 
               1        0.6     0.8];

%% load corr barcharts avg data for each channel & convert to network graphs
%{
% set up network nodes' coordinates, X and Y; # of nodes = num_band;
phi = (0:num_band-1).*2*pi/num_band + 0.5*pi;
X = cos(phi); Y = sin(phi);

% specify the prefix of data file names
data_prefix = 'barchart_cellarray_allsub_alltest_';

% loop through all channels
for chn_idx = 1:num_channel
    data_filename = strcat(data_prefix, chn_labels{chn_idx},'.mat');
    load([input_path, data_filename], 'corr_ratios_grp');
                       
    figure(chn_idx)
    
    % loop through all physiological stages
    for s_idx=1:num_stage  

        % extract group avg correlations (negative & positive)
        corr_ratio_avg = corr_ratios_grp(1:2,:,s_idx); 
        corr_ratio_err = corr_ratios_grp(3:4,:,s_idx); 

        % construct correlation networks (1=negative corr; 2=positive corr)
        for FLAG_corr = 1:2
            link_width = 8*(abs(corr_ratio_avg(FLAG_corr,:)) + 1e-6);  % adding eps=1e-6 to avoid zero link width
            link_color = 1 - abs(corr_ratio_avg(FLAG_corr,:));

            subplot(2,num_stage,s_idx+(2-FLAG_corr)*num_stage);     % row1=positive corr; row2=negative corr
            hold on;

            p_idx = 0;  % freq band index/counter
            for i = 1:num_band
                for j = i+1:num_band
                    p_idx = p_idx+1;
                    % plot network graph (xcorr probs as links; freq bands as nodes) 
                    plot([X(i) X(j)], [Y(i) Y(j)], ...
                        'linewidth',link_width(p_idx), ...
                        'color',link_color(p_idx).*ones(1,3))
                end
            end

            % plot nodes labels
            for i=1:num_band
                plot(X(i),Y(i), ...
                    'MarkerSize',14,'Marker','o','MarkerEdgeColor','k','linewidth',2);
                plot(X(i),Y(i),...
                    'MarkerSize',45,'Marker','.','MarkerEdgeColor',markercolor(i,:));

                text(X(i)-0.1,Y(i)+0.05,freqbands{i},'FontSize',12);

        %         plot(X(i),Y(i), ...
        %             'MarkerSize',75,'Marker','o','MarkerEdgeColor','k','linewidth',2);
        %         plot(X(i),Y(i),...
        %             'MarkerSize',250,'Marker','.','MarkerEdgeColor',markercolor(i,:));
        %         
        %         text(X(i)-0.1,Y(i)+0.05,label{i},'FontSize',45);
            end
            axis equal; 
            axis off; 

        end
    end
end
%}

%% Re-organize figure layouts
% set link color FLAG
link_color_FLAG = 0;    % 0=gray-scale; 1=stage-colormap

% set up network nodes' coordinates, X and Y; # of nodes = num_band;
phi = (0:num_band-1).*2*pi/num_band + 0.5*pi;
X = cos(phi); Y = sin(phi);
% X =1.5*X; Y = 1.5*Y;  % increase the diameter of the graphs

% specify the prefix of data file names
data_prefix = 'barchart_cellarray_allsub_alltest_';

% construct correlation networks (1=negative corr; 2=positive corr)
for FLAG_corr = 1:2
    
    % fig_handle = figure('Position', get(0, 'Screensize'));
    fig_handle = figure('Position', [100 100 1200 250]);
    % equally divide a figure window into num_stage=7 panels
    panelplots = panel();
    panelplots.pack('h', num_stage)
    % further divide each panel into 3-by-2 subpanels
    for s_idx=1:num_stage
        panelplots(s_idx).pack(3,2)
    end
    % adjust distances between panels
    panelplots.de.margin = 6;
%     panelplots.de.margin = 5;

    panelplots(1).marginleft = 0;
    for s_idx=2:num_stage
%         panelplots(s_idx).marginleft = 18;
        panelplots(s_idx).marginleft = 14;

    end
    
    % loop through all physiological stages
    for s_idx=1:num_stage  
%         figure(s_idx+(FLAG_corr-1)*num_stage)

        % loop through all channels
        for chn_idx = 1:num_channel
            data_filename = strcat(data_prefix, chn_labels{chn_idx},'.mat');
            load([input_path, data_filename], 'corr_ratios_grp');

            % extract group avg correlations (negative & positive)
            corr_ratio_avg = corr_ratios_grp(1:2,:,s_idx); 
            corr_ratio_err = corr_ratios_grp(3:4,:,s_idx); 

            link_width = 6*(abs(corr_ratio_avg(FLAG_corr,:)) + 1e-6);  % link thickness proportional to xcorr degree; 
                                                                       % adding eps=1e-6 to avoid zero link width
            link_color = 1 - abs(corr_ratio_avg(FLAG_corr,:));         % gray-scale link color proportional to xcorr degree
            
            % locate coordinates of networks at all channels
            row = ceil(chn_idx/2);  col = 2-mod(chn_idx,2);
            panelplots(s_idx,row,col).select();
            hold on;

            p_idx = 0;  % freq band index/counter
            for i = 1:num_band
                for j = i+1:num_band
                    p_idx = p_idx+1;
                    % plot network graph (xcorr probs as links; freq bands as nodes)                   
                    switch link_color_FLAG
                        case 0      % gray-scale 
                            plot([X(i) X(j)], [Y(i) Y(j)], ...
                                'linewidth',link_width(p_idx), ...
                                'color',link_color(p_idx).*ones(1,3))
                        case 1      % same colormap as histogram (colorcoded freq-pairs)                  
                            plot([X(i) X(j)], [Y(i) Y(j)], ...
                                'linewidth',link_width(p_idx), ...
                                'color',histcolor(p_idx,:))
                    end
                end
            end

            % plot nodes labels
            for i=1:num_band
                % use same colormap for each stage to color-code freq nodes
%                 plot(X(i),Y(i),'MarkerSize',50,'Marker','.','MarkerEdgeColor',labelcolor(s_idx,:));
                plot(X(i),Y(i),'MarkerSize',55,'Marker','.','MarkerEdgeColor',labelcolor(s_idx,:));

                % add greek letters for each freq node
%                 text(X(i)-0.1,Y(i)+0.05,freqbands{i},'FontSize',14);
                text(X(i)-0.175,Y(i)+0.05,freqbands{i},'FontSize',18);

                
                % use markercolor to color-code different freq nodes
                % plot(X(i),Y(i),'MarkerSize',50,'Marker','.','MarkerEdgeColor',markercolor(i,:));

                % Larger size plots
        %         plot(X(i),Y(i),'MarkerSize',250,'Marker','.','MarkerEdgeColor',markercolor(i,:));  
        %         text(X(i)-0.1,Y(i)+0.05,label{i},'FontSize',45);
        %         plot(X(i),Y(i),'MarkerSize',75,'Marker','o','MarkerEdgeColor','k','linewidth',2);
            end
            axis square; axis off; 
        end
        
    end
        
    % Save figures
    switch FLAG_corr
        case 1
            fig_name = 'network_graphs_negative';
        case 2
            fig_name = 'network_graphs_positive';
    end
%     savefig(fig_handle, [output_path, fig_name]);
%     saveas(fig_handle,[output_path, fig_name],'epsc')
%     exportgraphics(fig_handle,[output_path, fig_name,'.eps'],'Resolution',600)
    
end


%% calculate network weights (degree), link thickness & link color matrices

% and save them into excel sheets
for s_idx = 1: num_stage
    
    for chn_idx = 1:num_channel
        % define empty arrays for positive and negative correlation degrees
        pos_degree_mat = eye(num_band);
        neg_degree_mat = zeros(num_band);
        
        % load pairwise correlation for a given channel 
        data_filename = strcat(data_prefix, chn_labels{chn_idx},'.mat');
        load([input_path, data_filename], 'corr_ratios_grp');

        % allocate positive and negative degrees into matrices
        p_idx = 0;
        for row = 1:num_band
            for col = row+1:num_band
                p_idx = p_idx + 1;
                neg_degree_mat(row,col) = corr_ratios_grp(1,p_idx,s_idx);
                pos_degree_mat(row,col) = corr_ratios_grp(2,p_idx,s_idx);
            end
        end
        
        % symmetrify the matrices
        for col = 1:num_band
            for row = col+1:num_band
                neg_degree_mat(row,col) = neg_degree_mat(col,row);
                pos_degree_mat(row,col) = pos_degree_mat(col,row);
            end
        end
        
        % compute links thickness and link colors given degree matrices
        pos_link_thickness_mat = 6*(abs(pos_degree_mat) + 1e-6);
        neg_link_thickness_mat = 6*(abs(neg_degree_mat) + 1e-6);
        
        pos_link_color_mat = 1 - abs(pos_degree_mat);
        neg_link_color_mat = 1 - abs(neg_degree_mat);
        
        % convert these matrices into tables
        pos_degree_table = array2table(pos_degree_mat, "VariableNames",freqbands,"RowNames",freqbands); 
        neg_degree_table = array2table(neg_degree_mat, "VariableNames",freqbands,"RowNames",freqbands); 
        
        pos_link_thickness_table = array2table(pos_link_thickness_mat, "VariableNames",freqbands,"RowNames",freqbands); 
        neg_link_thickness_table = array2table(neg_link_thickness_mat, "VariableNames",freqbands,"RowNames",freqbands); 
        
        pos_link_color_table = array2table(pos_link_color_mat, "VariableNames",freqbands,"RowNames",freqbands); 
        neg_link_color_table = array2table(neg_link_color_mat, "VariableNames",freqbands,"RowNames",freqbands); 
        
        % save tables as excel sheets
        exel_loc = '../graph_data_excel_sheets/';
        excel_name = 'fig_5b.xlsx';
        file_name = strcat(exel_loc,excel_name);
        
        writetable(pos_degree_table, file_name, ... 
                   'Sheet',strcat(stages{s_idx},'_',chn_labels{chn_idx},'_pos_degree'),'WriteRowNames',true)
        writetable(neg_degree_table, file_name, ... 
                   'Sheet',strcat(stages{s_idx},'_',chn_labels{chn_idx},'_neg_degree'),'WriteRowNames',true)
        
        writetable(pos_link_thickness_table, file_name, ... 
                   'Sheet',strcat(stages{s_idx},'_',chn_labels{chn_idx},'_pos_link_thickness'),'WriteRowNames',true)
        writetable(neg_link_thickness_table, file_name, ... 
                   'Sheet',strcat(stages{s_idx},'_',chn_labels{chn_idx},'_neg_link_thickness'),'WriteRowNames',true)
        
        writetable(pos_link_color_table, file_name, ... 
                   'Sheet',strcat(stages{s_idx},'_',chn_labels{chn_idx},'_pos_link_color'),'WriteRowNames',true)
        writetable(neg_link_color_table, file_name, ... 
                   'Sheet',strcat(stages{s_idx},'_',chn_labels{chn_idx},'_neg_link_color'),'WriteRowNames',true)
    end
end
