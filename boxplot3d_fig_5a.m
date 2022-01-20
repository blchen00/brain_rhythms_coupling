% Generate network graphs from barcharts avg correlation profiles 
% Inputs:   input_path,  output_path
%           barcharts data: corr_ratios_grp, corr_ratios;
% 
% Outputs:  networks graphs 
% 
% 1/12/21: Bolun Chen
%
clear; clc; close all
%% Basic Settings & Parameters
% set input & output paths
input_path = '../data/barcharts_data/';       % where input .mat locates
output_path = '../raw_figs/';     % where to save figures
fig_name = 'xcorr_bars_evolutions';         % figure name

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
statelabel_short = {'R1';'WU';'EX';'CD';'CT1';'R2';'CT2'};
statelabel = {'Rest';'Warm Up';'Exercise';'Cool Down';'Cognitive Task';'Rest';'Cognitive Task'};

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

%% 3D bar chart plots showing evolution of correlation degree across physiologic states
% define an array for degrees of pairwise correlations across physiologic stages  
corr_degree = zeros(num_pair,num_stage);    

% specify the prefix of data file names
data_prefix = 'barchart_cellarray_allsub_alltest_';

chn_idx = 3;  % C3 channel
data_filename = strcat(data_prefix, chn_labels{chn_idx},'.mat');
load([input_path, data_filename], 'corr_ratios_grp');

% plot 3d bar charts with heights being the correlational degree
barwidth = 0.25;    % width of a 3d box
% x, y coordinates of region of box plots
xloc_min = 1; xloc_max = 3; xloc_separate = 4;  % separatrix between negative & positive boxes

barwidth = 0.6;    % width of a 3d box
xloc_min = 1; xloc_max = 7.25; xloc_separate = 9;  % separatrix between negative & positive boxes

xloc_min2 = xloc_max + 2*(xloc_separate-xloc_max);
xloc_max2 = xloc_min2 + (xloc_max-xloc_min);

yloc_min = 1; yloc_max = 15;
% yloc_separates = [5.5, 12.5];       % positive | mixed | negative correlations 
yloc_separates = [8.5, 20.5];       % positive | mixed | negative correlations 

% xloc = linspace(xloc_min,xloc_max,num_stage);
xloc1 = linspace(xloc_min,xloc_max,num_stage);
xloc2 = linspace(xloc_min2,xloc_max2,num_stage);
xloc = [xloc1, xloc2];

% yloc = linspace(yloc_min,yloc_max,num_pair);      % small intra-freq gap
yloc = [1:1.5:7,10:1.5:19,22:1.5:25];               % large intra-freq gap

[X,Y]= meshgrid(xloc,yloc);


figure('Position', [610 310 960 700])
hold on

% for negative and postive correlations (1=negative corr; 2=positive corr)
for FLAG_corr = 1:2

    % loop through all physiological stages to fill the array corr_degree
    for s_idx=1:num_stage  

        % extract group avg correlations (negative & positive)
        corr_ratio_avg = corr_ratios_grp(1:2,:,s_idx); 
        corr_ratio_err = corr_ratios_grp(3:4,:,s_idx); 

        corr_degree(:,s_idx) = corr_ratio_avg(FLAG_corr,:)';
    end
    
% plot negative & positive correlation bars
Z = (-1)^FLAG_corr.*corr_degree;        % degree of correlation is equal to 
                                        % the probability of obtaining 
                                        % certain value of xcorrelation
                                        % thus must be non-negative
% 3D bar plots (negative xcorr first, then postive xcorr)                                        
scatterbar3(X(:,1+(FLAG_corr-1)*num_stage:num_stage+(FLAG_corr-1)*num_stage),Y,Z,barwidth,labelcolor,histcolor)

end

% separating negative & positive correlation bars
fplot3(@(u) xloc_separate, @(u) u, @(u) 0, [min(yloc)-0.5 max(yloc)+1.5], 'LineStyle','--','LineWidth',1,'color','k')
fplot3(@(u) xloc_separate, @(u) max(yloc)+1.5, @(u) u, [0 1], 'LineStyle','--','LineWidth',1,'color','k')

% separating three categories of freq-pairs by solid lines
fplot3(@(u) u, @(u) yloc_separates(1), @(u) 0, [xloc_min-0.5 xloc_max2+0.5], 'LineStyle','-','LineWidth',1,'color','k')
fplot3(@(u) u, @(u) yloc_separates(2), @(u) 0, [xloc_min-0.5 xloc_max2+0.5], 'LineStyle','-','LineWidth',1,'color','k')
fplot3(@(u) xloc_max2+0.5, @(u) yloc_separates(1), @(u) u, [0 1], 'LineStyle','-','LineWidth',1,'color','k')
fplot3(@(u) xloc_max2+0.5, @(u) yloc_separates(2), @(u) u, [0 1], 'LineStyle','-','LineWidth',1,'color','k')

% adding texts to the figure
%{
text(0.5*(xloc_min+xloc_max)-0.9,max(yloc)+1.5,1.5, {'Degree of','anti-correlation'},...
    'HorizontalAlignment','center','FontSize',18,'FontWeight','bold','Rotation',40)
text(0.5*(xloc_min2+xloc_max2)-2.2,max(yloc)+1.5,1.5, {'Degree of','pos-correlation'}, ...
    'HorizontalAlignment','center','FontSize',18,'FontWeight','bold','Rotation',40)

% add frames around texts
patch('Position',[0,max(yloc)+1.5,1,0.5],'FaceColor','y','EdgeColor','k','LineWidth',2)


text(xloc_max2+0.5,3.5,1.5, {'Anti-correlated','brain rhythms'},...
    'HorizontalAlignment','center','FontSize',18,'FontWeight','bold','Rotation',-21)
text(xloc_max2+0.5,13.25,1.5, {'Mix-correlated','brain rhythms'},...
    'HorizontalAlignment','center','FontSize',18,'FontWeight','bold','Rotation',-21)
text(xloc_max2+0.5,22.6,1.6, {'Positive-','correlated','brain rhythms'},...
    'HorizontalAlignment','center','FontSize',18,'FontWeight','bold','Rotation',-21)
%}


% Color pairlabel (YTickLabel) with histcolormap
for i = 1:num_pair
    pairlabel{i} = sprintf('\\color[rgb]{%f,%f,%f}%s', histcolor(i,:), pairlabel{i});
end

% Color stagebale (XTickLabel) with labelcolormap
for i = 1:num_stage
%     statelabel{i} = sprintf('\\color[rgb]{%f,%f,%f}%s', labelcolor(i,:), statelabel{i});
        statelabel_short{i} = sprintf('\\color[rgb]{%f,%f,%f}%s', labelcolor(i,:), statelabel_short{i});

end

% figure detail setup
% axis equal
pbaspect([10.6 18.3 2])
pbaspect([20 32 4])

view(-30,30)
view(-55,35)
%{
box on; 
set(gca,'XGrid','off','YGrid','off','ZGrid','on', ...
    'XTick',xloc,'XTickLabel',[statelabel statelabel],...
    'YTick',yloc,'YTickLabel',pairlabel,...
    'XLim',[xloc_min-0.5 xloc_max2+0.5],...
    'YLim',[yloc_min-0.5 max(yloc)+1.5],...% 'YLim',[yloc_min-0.5 yloc_max+0.5]
    'ZLim',[0 0.9],...
    'FontSize',16.5,'FontWeight','bold',...
    'XTickLabelRotation',-40,'YTickLabelRotation',16);
%}

box on; 
set(gca,'XGrid','off','YGrid','off','ZGrid','on', ...
    'XTick',xloc,'XTickLabel',[statelabel_short statelabel_short],...
    'YTick',yloc,'YTickLabel',pairlabel,...
    'XLim',[xloc_min-0.5 xloc_max2+0.5],...
    'YLim',[yloc_min-0.5 max(yloc)+1.5],...% 'YLim',[yloc_min-0.5 yloc_max+0.5]
    'ZLim',[0 0.9],...    
    'FontSize',20,'FontWeight','bold',...
    'XTickLabelRotation',-24,'YTickLabelRotation',44);
ax = gca; ax.YAxis.FontSize = 28; ax.XAxis.FontSize = 18;


% Save figures

% savefig(gcf, [output_path, fig_name]);
% saveas(gcf,[output_path, fig_name],'epsc')

%% save bar data into excel sheets
% both fig_4 and fig_5a use the same data arrays
neg_corr_array = squeeze(corr_ratios_grp(1,:,:))*(-1);
pos_corr_array = squeeze(corr_ratios_grp(2,:,:));

pairlabel = cell(1,num_pair);
count = 0;
for i=1:num_band-1
    for j=i+1:num_band
        pair = strcat(freqbands{i},  '-', freqbands{j});
        count = count+1;
        pairlabel{count} = pair;
    end
end      

% convert arrays to tables
neg_corr_table = array2table(neg_corr_array, ...
                             'VariableNames', {stages{:}}, ...
                             'RowNames',{pairlabel{:}}); 
pos_corr_table = array2table(pos_corr_array, ...
                             'VariableNames', {stages{:}}, ...
                             'RowNames',{pairlabel{:}}); 

% save these tables as Excel spreadsheets for fig_4
exel_loc = '../graph_data_excel_sheets/';
excel_name = 'fig_4a.xlsx';
file_name = strcat(exel_loc,excel_name);

writetable(pos_corr_table, file_name, 'Sheet','pos_degree','WriteRowNames',true)
writetable(neg_corr_table, file_name, 'Sheet','neg_degree','WriteRowNames',true)

% save these tables as Excel spreadsheets for fig_5a
exel_loc = '../graph_data_excel_sheets/';
excel_name = 'fig_5a.xlsx';
file_name = strcat(exel_loc,excel_name);

writetable(pos_corr_table, file_name, 'Sheet','pos_degree','WriteRowNames',true)
writetable(neg_corr_table, file_name, 'Sheet','neg_degree','WriteRowNames',true)


%% 
% generate the raw version of fig_4a
figure
colormap jet
% colormap turbo
% colormap parula
subplot(211)
imagesc(pos_corr_array)
colorbar
caxis([0 1])
xlabel('time')
title('degree of positive correlation')

subplot(212)
imagesc(neg_corr_array)
colorbar
caxis([0 1])
xlabel('time')
title('degree of negative correlation')

% generate the version used in fig_4a
figure
colormap jet
% colormap turbo
subplot(211)
imagesc(pos_corr_array)
colorbar
caxis([0 1])
set(gca,'XTickLabel',[],'XTick',[],'YTick',[])
% xlabel('time')
% title('degree of positive correlation')

subplot(212)
imagesc(neg_corr_array)
colorbar
caxis([0 1])
set(gca,'XTickLabel',[],'XTick',[],'YTick',[])
% xlabel('time')
% title('degree of negative correlation')




