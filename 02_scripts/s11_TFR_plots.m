clc; clear all; close all;
%% prepare
addpath('C:\Users\chris\Desktop\EPOC\EPOC-Kiel\04_software\fieldtrip-20250106');

load("C:\Users\chris\Desktop\EPOC\EPOC-Kiel\01_data\topoplot_workspace.mat"); 

% path to folder with ICA weights
path = "C:\Users\chris\Desktop\EPOC\EPOC-Kiel\01_data\01_individual_data_after_ICA";
% list of all .mat files
data = dir(fullfile(path, '*.mat'));

%% build groups
T_merged        = readtable("C:\Users\chris\Desktop\EPOC\EPOC-Kiel\01_data\00_bids\participants.tsv", 'Delimiter', '\t', 'FileType','text');


used_ids        = cellfun(@(x) x(1:end-4), {data.name}, 'UniformOutput', false);

valid_idx       = ~cellfun(@isempty, TFR);
TFR             = TFR(valid_idx);
used_ids        = used_ids(valid_idx);

TFR([49])       = [];
used_ids([49])  = [];

% Mapping for T_merged
[~, idx]        = ismember(used_ids, T_merged.participant_id);
group_labels    = T_merged.group(idx);

idx_with        = strcmp(group_labels, 'withPCS');
idx_without     = strcmp(group_labels, 'withoutPCS');

TFR_with        = TFR(idx_with);
TFR_without     = TFR(idx_without);

%% Grand Average per group 
GA_with         = ft_freqgrandaverage([], TFR_with{:});
GA_without      = ft_freqgrandaverage([], TFR_without{:});

%% settings 
freqband        = [2 15];                                                   % Hz
timewin         = [-1.0 1.0];                                               % seconds
roi             = 'all';

% colors
zlimBoth = get_common_zlim_both(GA_with, GA_without, freqband, timewin, roi);

%% ===== Plot Figure =====
figure('Color','w');
set(gcf,'Units','pixels','Position',[100 100 2500 2000]);
set(gcf,'PaperPositionMode','auto');

tlo = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

% colors
pink = [204 121 167] / 255;   % #CC79A7
wine = [130 50 94]   / 255;   % #82325E

% --- group WITH ---
nexttile;
plot_single_TFR_subplot(GA_with, ...
    'Self-Reported Cognitive Symptoms', freqband, timewin, roi, zlimBoth);
title('Self-Reported Cognitive Symptoms', ...
    'FontSize', 40, 'FontWeight','bold', 'Color', pink);

% --- group WITHOUT ---
nexttile;
plot_single_TFR_subplot(GA_without, ...
    'No Self-Reported Cognitive Symptoms', freqband, timewin, roi, zlimBoth);
title('No Self-Reported Cognitive Symptoms', ...
    'FontSize', 40, 'FontWeight','bold', 'Color', wine);

annotation('textbox', [0.05 0.965 0.6 0.03], ...
    'String', 'Neumann - Neurophysiological Correlates of Cognitive Symptoms in Post−COVID−19', ...
    'HorizontalAlignment', 'left', ...
    'EdgeColor', 'none', ...
    'FontSize', 11, ...
    'FontWeight', 'normal')

exportgraphics(gcf, 'C:\Users\chris\Desktop\EPOC\EPOC-Kiel\03_figures\07_eeg_peaks\TFR.tif', 'Resolution', 500);

%% ===== functions =====

function plot_single_TFR_subplot(GA, titleText, freqband, timewin, roi, zlimFix)

    hasChan = isfield(GA,'dimord') && contains(GA.dimord,'chan');
    if ischar(roi) && strcmp(roi,'all')
        channels = 'all';
    else
        channels = roi;
    end

    cfg = [];
    cfg.parameter   = 'powspctrm';
    cfg.xlim        = timewin;
    cfg.ylim        = freqband;
    cfg.zlim        = zlimFix;
    cfg.interactive = 'no';
    cfg.figure      = 'gca';

    if hasChan
        cfg.channel     = channels;
        cfg.avgoverchan = 'yes';
    end

    ft_singleplotTFR(cfg, GA);

    hold on

    % stimulus line
    xline(0,'k','LineWidth',2);

    % colors
    thetaColor = [230 159 0] / 255;   % #E69F00
    alphaColor = [192 133 0] / 255;   % #C08500

    % theta 4–7 Hz
    yline(4,'Color',thetaColor,'LineWidth',4);
    yline(7,'Color',thetaColor,'LineWidth',4);

    % tlpha 8–13 Hz
    yline(8, 'Color',alphaColor,'LineWidth',4);
    yline(13,'Color',alphaColor,'LineWidth',4);

    hold off

    % labels bigger
    xlabel('Time (s)','FontSize',22);
    ylabel('Frequency (Hz)','FontSize',22);

    % text bigger
    ax = gca;
    ax.FontSize = 22;
end


function zlimCommon = get_common_zlim_both(GA1, GA2, freqband, timewin, roi)

    [zmins, zmaxs] = deal([]);
    GAall = {GA1, GA2};

    for k = 1:2
        G = GAall{k};
        hasChan = isfield(G,'dimord') && contains(G.dimord,'chan');

        cfgs = [];
        cfgs.frequency = freqband;
        cfgs.latency   = timewin;

        if hasChan
            if ischar(roi) && strcmp(roi,'all')
                cfgs.channel = 'all';
            else
                cfgs.channel = roi;
            end
            cfgs.avgoverchan = 'yes';
        end

        sub = ft_selectdata(cfgs, G);
        if isfield(sub,'powspctrm') && ~isempty(sub.powspctrm)
            zmins(end+1) = nanmin(sub.powspctrm(:));
            zmaxs(end+1) = nanmax(sub.powspctrm(:));
        end
    end

    zlimCommon = [min(zmins) max(zmaxs)];
end
