clc; clear all; close all;

%% prepare everything

% path to folder with ICA weights
path                = "C:\Users\chris\Desktop\EPOC\EPOC-Kiel\01_data\01_individual_data_after_ICA";

% list of all .mat files
data = dir(fullfile(path, '*.mat'));

MAIN                = 'C:\Users\chris\Desktop\EPOC\EPOC-Kiel\';                       % paste path to EPOC folder
cd(MAIN);                                                                   % go into EPOC folder

epocPath            = genpath(MAIN);
fieldtripPath       = genpath('C:\Users\chris\Desktop\EPOC\EPOC-Kiel\04_software\fieldtrip-20250106');
epocPath            = strrep(epocPath, fieldtripPath, '');                  % add fieldtrip seperately from the rest 

addpath(epocPath);                                                          % add EPOC 
addpath('C:\Users\chris\Desktop\EPOC\EPOC-Kiel\04_software\fieldtrip-20250106');      % and fieldtrip 

PATHIN      = [MAIN '01_data\00_bids\'];
PATHOUT     = [MAIN '01_data\'];


% prepare patients 
patient             = dir(PATHIN);
patient([1,2],:)    = [];

isVP                = startsWith({patient.name}, "sub");
patient             = patient(isVP);

% prepare layout
cfg = [];
load('BC-128-pass-lay.mat');
cfg.elec            = elec;
layout              = ft_prepare_layout(cfg);


% Match the updated main EEG cohort before running the time-consuming TFR analysis.
main_results = load([PATHOUT '02_fooof_pvt' filesep 'fooofed.mat'], ...
    'results_pre', 'results_post');
main_ids = string({main_results.results_pre.id})';
main_ids = main_ids(strlength(main_ids) > 0);
post_ids = string({main_results.results_post.id})';
post_ids = post_ids(strlength(post_ids) > 0);
assert(isequal(sort(main_ids), sort(post_ids)), ...
    'Pre/Post participant IDs in fooofed.mat do not match.');
assert(numel(unique(main_ids)) == numel(main_ids), ...
    'Duplicate IDs in fooofed.mat.');

% Set to false when the TFR is already saved and only the plot should be recreated.
recomputeTFR = true;
workspace_file = [PATHOUT 'topoplot_workspace.mat'];

if recomputeTFR

TFR = cell(numel(data), 1);

%% Loop over every participant

for i               = 1:length(data)
    %% read raw data 

    datapath        = fullfile(path, data(i).name);
    load(datapath);
    sub             = data(i).name(1:end-4);

    % Reproduce the participant set of the updated peak analysis.
    if ~ismember(string(sub), main_ids)
        continue;
    end

    %% delete bad channels
    idx             = find(strcmp({patient.name}, sub));                    % find the Index in patient
    cfg             = [];
    cfg.channel     = setdiff(data_r.label, bad_channels_all{idx});         % remove noisy channels
    data_c          = ft_selectdata(cfg, data_r);

    %% apply filters

    cfg = [];
    cfg.demean      = 'yes';                                                % remove DC offset
    cfg.hpfilter    = 'yes';
    cfg.hpfreq      = 1;                                                    % high-pass filter, cutting everything under 1 Hz
    cfg.hpfilttype  = 'firws';
    cfg.lpfilter    = 'yes';
    cfg.lpfreq      = 90;                                                   % low-pass filter, cutting everything over 90 Hz 
    cfg.lpfilttype  = 'firws';
    cfg.bsfilter    = 'yes'; 
    cfg.bsfreq      = [49 51];                                              % remove 50 Hz frequencies from outlet
    cfg.bsfilttype  = 'firws';
    cfg.pad         = 'nextpow2';
    data_p          = ft_preprocessing(cfg, data_c);

    %% segment data

    % search for pvt data
    cd([PATHIN sub filesep 'eeg']);
    indat                   = dir(fullfile(PATHIN, patient(idx).name, 'eeg', '*pvt_eeg.vhdr'));   

    % search for events
    cfg                     = [];
    cfg.dataset             = [PATHIN sub filesep 'eeg' filesep indat.name];
    events                  = ft_read_event(cfg.dataset);

    % prepare ft_definetrial
    cfg.trialdef.eventtype  = 'Markers';
    cfg.trialdef.eventvalue = 'Trial start';
    cfg.trialdef.prestim    = 1.2;
    cfg.trialdef.poststim   = 1.2;
    
    try
        cfg                 = ft_definetrial(cfg); 
    catch ME
        % skip loop 
        warning('Skipping subject %d due to error: %s', i, ME.message);
        continue; 
    end

    %% Remove practice trials and invalid reaction times

    % Get trial start and button press markers
    values = string({events.value});

    start_samples = [events(values == "Trial start").sample]';
    bp_samples    = [events(values == "BP").sample]';

    % Match EEG trials to their original stimulus markers
    trl_onsets = cfg.trl{:,1} - cfg.trl{:,3};
    nTrials = size(cfg.trl, 1);

    marker_idx = zeros(nTrials, 1);

    for t = 1:nTrials
        matches = find(abs(start_samples - trl_onsets(t)) <= 1);

        assert(numel(matches) == 1, ...
            'No unique marker match for trial %d in %s', t, sub);

        marker_idx(t) = matches;
    end

    % Check that each marker is used only once
    assert(numel(unique(marker_idx)) == nTrials, ...
        'Duplicate marker assignment for %s', sub);

    % Calculate reaction times
    rt = NaN(nTrials, 1);

    for t = 1:nTrials

        m = marker_idx(t);

        % Find next stimulus marker
        if m < numel(start_samples)
            next_start = start_samples(m+1);
        else
            next_start = Inf;
        end

        % Find first button press after stimulus
        bp = bp_samples(bp_samples > start_samples(m) & ...
                        bp_samples < next_start);

        if ~isempty(bp)
            rt(t) = (bp(1) - start_samples(m)) / data_p.fsample;
        end
    end

    % Identify valid trials using ORIGINAL trial numbers
    valid_trials = marker_idx > 8 & ...
                   rt >= 0.1 & rt <= 1.0;

    % Display trial counts
    fprintf('%s:\n', sub);
    fprintf('  Stimulus markers: %d\n', numel(start_samples));
    fprintf('  EEG trials before exclusion: %d\n', nTrials);
    fprintf('  Practice trials excluded: %d\n', sum(marker_idx <= 8));
    fprintf('  RT < 100 ms: %d\n', sum(rt(marker_idx > 8) < 0.1));
    fprintf('  RT > 1000 ms: %d\n', sum(rt(marker_idx > 8) > 1.0));
    fprintf('  Missing button presses: %d\n', sum(isnan(rt(marker_idx > 8))));
    fprintf('  Remaining trials: %d\n', sum(valid_trials));

    % Remove excluded trials
    cfg.trl = cfg.trl(valid_trials, :);

    % Segment remaining trials
    data_t = ft_redefinetrial(cfg, data_p);

    %% delete bad components
    
    cfg             = [];
    cfg.component   = rejected_comps(:)'; 
    comp.label      = arrayfun(@(x) sprintf('IC%d', x), 1:size(comp.unmixing,1), 'UniformOutput', false);

    data_ti         = ft_rejectcomponent(cfg,comp, data_t); 

    %% remove eog

    cfg             = [];
    cfg.channel     = 'all';
    cfg.channel     = ft_channelselection({'all', '-31', '-32'}, data_ti.label);
    data_ti         = ft_selectdata(cfg, data_ti);

    %% remove bad epochs

    orig_labels             = data_ti.label;
    
    % convert to eeglab
    EEG                     = fieldtrip2eeglab(data_ti);
    
    % detect and delete artefacts with jointprob
    % Check for NaNs or Infs
    if any(isnan(EEG.data(:))) || any(isinf(EEG.data(:)))
        warning('data contains NaNs or Infs – pop_jointprob skipped.');
        continue;
    else
        [EEG, bad_epochs]   = pop_jointprob(EEG, 1, 1:EEG.nbchan, 3, 3, 0, 1, 0);
    end

    
    % convert back to fieldtrip
    data_tc                 = eeglab2fieldtrip(EEG, 'preprocessing');
    data_tc.label           = orig_labels;

    %% check if there are trials left

    if isempty(data_tc.trial)
        fprintf('no trials left for particpant %s. Skipping...\n', sub); 
        no_trials_left{i} = sub;
        continue; 
    end

    %% Repair deleted channels

    % prepare Neighbors
    cfg                     = []; 
    cfg.method              = 'distance';                                   % how should the neighbors be selected?
    cfg.neighbourdist       = 0.5; 
    cfg.elec                = elec;
    neigh                   = ft_prepare_neighbours(cfg);                   % between 5 and 10 neighbors is a good value, always good to check!
        
    % Check Neighbors
%     cfg = [];
%     cfg.neighbours = neigh; % what neighbor-structure
%     cfg.elec = elec;
%     
%     ft_neighbourplot(cfg)

    cfg                     = [];
    cfg.missingchannel      = bad_channels_all{idx};
    cfg.method              = 'spline';
    cfg.neighbours          = neigh;
    cfg.elec                = elec;


    try
    data_n                  = ft_channelrepair(cfg, data_tc);
    catch
        disp('No missing channels to repair. Skipping ft_channelrepair.');
        data_n = data_tc; 
    end

    %% Re-Referencing
    
    cfg                 = [];
    cfg.reref           = 'yes';
    cfg.refchannel      = find(~strcmpi(data_n.label,'Gnd'));               % Take all channels
    cfg.refmethod       = 'avg';                                            % Take the average
    
    data_erp            = ft_preprocessing(cfg,data_n);


    %% Time-Frequency Analysis

    cfg = [];
    cfg.method          = 'mtmconvol';
    cfg.output          = 'pow';
    cfg.taper           = 'hanning';
    cfg.foi             = 2:0.5:15;
    cfg.t_ftimwin       = 2 ./ cfg.foi;
    cfg.toi             = -1.2:0.01:1.2;
    cfg.pad             = 'nextpow2';
    cfg.padtype         = 'zero';
    TFR{i}              = ft_freqanalysis(cfg, data_erp);

end

%% Save TFR data before plotting

file_ids = cellfun(@(x) x(1:end-4), {data.name}, 'UniformOutput', false);
valid_idx = ~cellfun(@isempty, TFR);
TFR = TFR(valid_idx);
used_ids = string(file_ids(valid_idx))';

% Include exactly the participants in the revised EEG main analysis.
keep_idx = ismember(used_ids, main_ids);
TFR = TFR(keep_idx);
used_ids = used_ids(keep_idx);

assert(numel(used_ids) == numel(main_ids) && ...
       numel(unique(used_ids)) == numel(used_ids), ...
       'The TFR cohort does not match fooofed.mat.');

% Keep the previous channel-wise spectra in this workspace intact.
TFR_ids = used_ids;
if isfile(workspace_file)
    save(workspace_file, 'TFR', 'TFR_ids', '-append');
else
    save(workspace_file, 'TFR', 'TFR_ids', '-v7.3');
end

fprintf('TFR saved for %d participants: %s\n', numel(TFR_ids), workspace_file);

else
    assert(isfile(workspace_file), 'TFR workspace file does not exist.');
    tfr_saved = load(workspace_file, 'TFR', 'TFR_ids');
    TFR = tfr_saved.TFR;
    used_ids = string(tfr_saved.TFR_ids);
end

% Confirm that the saved/recomputed TFR still matches the main EEG cohort.
assert(numel(used_ids) == numel(main_ids) && ...
       all(ismember(used_ids, main_ids)) && ...
       numel(unique(used_ids)) == numel(used_ids), ...
       'Saved TFR participant IDs do not match fooofed.mat.');

%% Build groups

T_merged = readtable([PATHIN 'participants.tsv'], ...
    'Delimiter', '\t', 'FileType', 'text', 'TextType', 'string');

[found, idx] = ismember(used_ids, string(T_merged.participant_id));
assert(all(found), 'Some TFR IDs were not found in participants.tsv.');

group_labels = lower(strtrim(string(T_merged.group_Fatigue(idx))));
idx_with = group_labels == "clinically relevant fatigue";
idx_without = group_labels == "no clinically relevant fatigue";

unknown = ~(idx_with | idx_without);
if any(unknown)
    warning('%d participants have unexpected or missing fatigue group labels and are not plotted.', ...
        sum(unknown));
    disp(table(used_ids(unknown), group_labels(unknown), ...
        'VariableNames', {'participant_id','group_Fatigue'}));
end

fprintf('Fatigue group: %d; no-fatigue group: %d; unclassified: %d\n', ...
    sum(idx_with), sum(idx_without), sum(unknown));
assert(any(idx_with) && any(idx_without), ...
    'At least one fatigue group is empty. Inspect group_Fatigue in participants.tsv.');

TFR_with = TFR(idx_with);
TFR_without = TFR(idx_without);

%% Grand Average per group

GA_with    = ft_freqgrandaverage([], TFR_with{:});
GA_without = ft_freqgrandaverage([], TFR_without{:});

%% Settings

freqband = [2 15];
timewin  = [-1.0 1.0];
roi      = 'all';

% Use the same color limits for both panels.
zlimBoth = get_common_zlim_both(GA_with, GA_without, freqband, timewin, roi);

%% Plot Figure

set(groot,'defaultAxesFontName','Arial')
set(groot,'defaultTextFontName','Arial')

figure('Color','w');
set(gcf,'Units','pixels','Position',[100 100 1800 1400]);
set(gcf,'PaperPositionMode','auto');

tlo = tiledlayout(2,1,...
    'TileSpacing','compact',...
    'Padding','compact');

% Colors
pink = [204 121 167] / 255;
wine = [130 50 94] / 255;

%% Group WITH

ax(1) = nexttile;

plot_single_TFR_subplot(GA_with, ...
    'Clinically significant fatigue', ...
    freqband, timewin, roi, zlimBoth);

title('Fatigue',...
    'FontName','Arial',...
    'FontSize',15,...
    'FontWeight','bold',...
    'Color',pink);

%% Group WITHOUT

ax(2) = nexttile;

plot_single_TFR_subplot(GA_without, ...
    'No clinically significant fatigue', ...
    freqband, timewin, roi, zlimBoth);

title('No Fatigue',...
    'FontName','Arial',...
    'FontSize',15,...
    'FontWeight','bold',...
    'Color',wine);

%% Formatting

for k = 1:2
    ax(k).FontName = 'Arial';
    ax(k).FontSize = 10;
    ax(k).XAxis.FontSize = 10;
    ax(k).YAxis.FontSize = 10;
    ax(k).XLabel.FontSize = 11;
    ax(k).YLabel.FontSize = 11;
    ax(k).LineWidth = 1;
end

%% Panel labels

letters = {'A','B'};

for k = 1:2
    annotation('textbox',...
        [ax(k).Position(1)-0.015,...
         ax(k).Position(2)+ax(k).Position(4)-0.005,...
         0.03,...
         0.03],...
        'String',letters{k},...
        'LineStyle','none',...
        'FontName','Arial',...
        'FontSize',15,...
        'FontWeight','bold',...
        'Margin',0);
end

%% Colorbar formatting

cb = findall(gcf,'Type','ColorBar');

for k = 1:length(cb)
    cb(k).FontName = 'Arial';
    cb(k).FontSize = 10;
    cb(k).Label.FontName = 'Arial';
    cb(k).Label.FontSize = 11;
end

%% Export

set(gcf, 'Units', 'centimeters', 'Position', [5 5 18 22])

exportgraphics(gcf,...
'C:\Users\chris\Desktop\EPOC\EPOC-Kiel\03_figures\06_eeg_peaks\TFR_2Hz_RTfiltered.png',...
'Resolution',600);

%% Functions

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

    % Stimulus line
    xline(0,'k','LineWidth',2);

    % Colors
    thetaColor = [230 159 0] / 255;   % #E69F00
    alphaColor = [192 133 0] / 255;   % #C08500

    % Theta 4-7 Hz
    yline(4,'Color',thetaColor,'LineWidth',4);
    yline(7,'Color',thetaColor,'LineWidth',4);

    % Alpha 8-12 Hz, consistent with the main analysis
    yline(8, 'Color',alphaColor,'LineWidth',4);
    yline(12,'Color',alphaColor,'LineWidth',4);

    hold off

    xlabel('Time (s)','FontSize',22);
    ylabel('Frequency (Hz)','FontSize',22);

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
