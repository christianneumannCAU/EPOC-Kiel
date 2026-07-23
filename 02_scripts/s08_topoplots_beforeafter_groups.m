clc; clear all; close all;

%% prepare everything

% path to folder with ICA weights
path                = "C:\Users\chris\Desktop\EPOC\EPOC-Kiel\01_data\01_individual_data_after_ICA";

% list of all .mat files
data = dir(fullfile(path, '*.mat'));

MAIN                = 'C:\Users\chris\Desktop\EPOC\EPOC-Kiel\';             % paste path to EPOC folder
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

corr_power_pre = cell(1,95);
corr_power_post = cell(1,95);

freq_pre = cell(95, 1);
freq_post = cell(95, 1);

%% Loop over every participant

for i               = 1:length(data)
    %% read raw data 

    datapath        = fullfile(path, data(i).name);
    load(datapath);
    sub             = data(i).name(1:end-4);

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

    data_t                  = ft_redefinetrial(cfg, data_p);

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
    % Check for NaNs oder Infs
    if any(isnan(EEG.data(:))) || any(isinf(EEG.data(:)))
        warning('data contains NaNs or Infs – pop_jointprob skipped.');
        continue;
    else
        [EEG, bad_epochs]   = pop_jointprob(EEG, 1, 1:EEG.nbchan, 3, 3, 0, 1, 0);
    end

    
    % convert back to fieldtrip
    data_tc                 = eeglab2fieldtrip(EEG, 'preprocessing');
    data_tc.label           = orig_labels;

    % How many trials per participant? 
    all_bad_trials{i}       = sum(EEG.reject.rejjp);

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


    %% ===== Time-Frequency Analysis =====
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

    %% average pre
    cfg             = [];
    cfg.latency     = [-1 0];
    data_pre        = ft_selectdata(cfg, data_erp);
    
    % average trials
    cfg             = [];
    data_timelock   = ft_timelockanalysis(cfg, data_pre);

    %% aperiodic pre 

    cfg                         = [];
    cfg.fooof.aperiodic_mode    = 'knee';
    cfg.foilim                  = [1 40]; 
    cfg.tapsmofrq               = 4; 
    cfg.pad                     = 4;
    cfg.method                  = 'mtmfft';
    cfg.output                  = 'fooof_aperiodic';

    fractal                     = ft_freqanalysis(cfg, data_timelock);

    cfg.output                  = 'pow';
    original                    = ft_freqanalysis(cfg, data_timelock);

    % subtract the fractal component from the power spectrum
    cfg                         = [];
    cfg.parameter               = 'powspctrm';
    cfg.operation               = 'x2-x1';
    oscillatory                 = ft_math(cfg, fractal, original);

    freq_pre{i} = oscillatory;

    %% average post

    cfg                         = [];
    cfg.latency = [0 1];
    data_post                   = ft_selectdata(cfg, data_erp);
    
    % average trials
    cfg                         = [];
    data_timelock               = ft_timelockanalysis(cfg, data_post);

    %% aperiodic post

    cfg                         = [];
    cfg.fooof.aperiodic_mode    = 'knee';
    cfg.foilim                  = [1 40]; 
    cfg.tapsmofrq               = 4; 
    cfg.pad                     = 4;
    cfg.method                  = 'mtmfft';
    cfg.output                  = 'fooof_aperiodic';

    fractal                     = ft_freqanalysis(cfg, data_timelock);

    cfg.output                  = 'pow';
    original                    = ft_freqanalysis(cfg, data_timelock);

    % subtract the fractal component from the power spectrum
    cfg                         = [];
    cfg.parameter               = 'powspctrm';
    cfg.operation               = 'x2-x1';
    oscillatory                 = ft_math(cfg, fractal, original);

    freq_post{i}                = oscillatory; 

end

save('C:\Users\chris\Desktop\EPOC\EPOC-Kiel\01_data\topoplot_workspace.mat');

%% ------------------------------------------------------------
%  THETA
% ------------------------------------------------------------

load("C:\Users\chris\Desktop\EPOC\EPOC-Kiel\01_data\topoplot_workspace.mat");
T_merged = readtable("C:\Users\chris\Desktop\EPOC\EPOC-Kiel\01_data\00_bids\participants.tsv", 'Delimiter', '\t', 'FileType','text');

used_ids = cellfun(@(x) x(1:end-4), {data.name}, 'UniformOutput',false);

valid_idx = ~cellfun(@isempty,freq_pre) & ~cellfun(@isempty,freq_post);
freq_pre  = freq_pre(valid_idx);
freq_post = freq_post(valid_idx);
used_ids  = used_ids(valid_idx);

freq_pre([49,60,70]) = [];
freq_post([49,60,70]) = [];
used_ids([49,60,70])  = [];

[~,idx] = ismember(used_ids,T_merged.participant_id);
group_labels = T_merged.group_Fatigue(idx);

idx_with    = strcmp(group_labels,'clinically relevant fatigue');
idx_without = strcmp(group_labels,'no clinically relevant fatigue');

freq_pre_with    = freq_pre(idx_with);
freq_pre_without = freq_pre(idx_without);
freq_post_with   = freq_post(idx_with);
freq_post_without= freq_post(idx_without);

GA_pre_with       = ft_freqgrandaverage([],freq_pre_with{:});
GA_pre_without    = ft_freqgrandaverage([],freq_pre_without{:});
GA_post_with      = ft_freqgrandaverage([],freq_post_with{:});
GA_post_without   = ft_freqgrandaverage([],freq_post_without{:});

col_pre   = [0 96 71]/255;      % #006047
col_post  = [0 158 115]/255;    % #009E73
col_self  = [204 121 167]/255;  % #CC79A7
col_no    = [130 50 94]/255;    % #82325E

col_theta_circle = [230 159 0]/255;  % #E69F00
col_alpha_circle = [192 133 0]/255;  % #C08500

cfg = [];
load('C:\Users\chris\Desktop\EPOC\EPOC-Kiel\BC-128-pass-lay.mat');
cfg.elec   = elec;
cfg.layout = ft_prepare_layout(cfg);
cfg.marker = 'on';
cfg.markersymbol = '.';
cfg.markersize = 3;
cfg.comment = 'no';
cfg.zlim = [-0.06 0.06];
cfg.xlim = [4 7];

set(groot,'defaultAxesFontName','Arial')
set(groot,'defaultTextFontName','Arial')

figure('Color','w','Position',[100 100 800 600]);
cfg.figure = gcf;
cfg.colormap = ft_colormap('balance',256);

ax(1) = subplot(2,2,1);
ft_topoplotER(cfg, GA_pre_with)
title({['\color[rgb]{' num2str(col_pre) '}Pre-Stimulus'], ...
       ['\color[rgb]{' num2str(col_self) '}Fatigue']}, ...
       'FontSize', 15, 'FontWeight', 'bold', 'Interpreter', 'tex')


ax(2) = subplot(2,2,2);
ft_topoplotER(cfg, GA_pre_without)
title({['\color[rgb]{' num2str(col_pre) '}Pre-Stimulus'], ...
       ['\color[rgb]{' num2str(col_no) '}No Fatigue']}, ...
       'FontSize', 15, 'FontWeight', 'bold', 'Interpreter', 'tex')


ax(3) = subplot(2,2,3);
ft_topoplotER(cfg, GA_post_with)
title({['\color[rgb]{' num2str(col_post) '}Post-Stimulus'], ...
       ['\color[rgb]{' num2str(col_self) '}Fatigue']}, ...
       'FontSize', 15, 'FontWeight', 'bold', 'Interpreter', 'tex')


ax(4) = subplot(2,2,4);
ft_topoplotER(cfg, GA_post_without)
title({['\color[rgb]{' num2str(col_post) '}Post-Stimulus'], ...
       ['\color[rgb]{' num2str(col_no) '}No Fatigue']}, ...
       'FontSize', 15, 'FontWeight', 'bold', 'Interpreter', 'tex')

% Adjust subplot positions
ax(1).Position = [0.08 0.56 0.28 0.30];
ax(2).Position = [0.34 0.56 0.28 0.30];
ax(3).Position = [0.08 0.16 0.28 0.30];
ax(4).Position = [0.34 0.16 0.28 0.30];

% Add panel labels
letters = {'A','B','C','D'};

for k = 1:4
    axes(ax(k))
    text(-0.05,1.05,letters{k},...
        'Units','normalized',...
        'FontName','Arial',...
        'FontWeight','bold',...
        'FontSize',16,...
        'HorizontalAlignment','left');
end

% Add colorbar
h = colorbar('Position',[0.70 0.16 0.018 0.74]);

ylabel(h,'Power (\muV^2)',...
    'FontName','Arial',...
    'FontSize',16,...
    'FontWeight','bold');

set(h,...
    'FontName','Arial',...
    'FontSize',14);

% Export
set(gcf,'Units','pixels','Position',[100 100 2500 1800])

exportgraphics(gcf,...
'C:\Users\chris\Desktop\EPOC\EPOC-Kiel\03_figures\06_eeg_peaks\topo_theta.pdf',...
'Resolution',600)

%% ------------------------------------------------------------
%  ALPHA
% ------------------------------------------------------------

cfg = [];
load('C:\Users\chris\Desktop\EPOC\EPOC-Kiel\BC-128-pass-lay.mat');
cfg.elec   = elec;
cfg.layout = ft_prepare_layout(cfg);
cfg.marker = 'on';
cfg.markersymbol = '.';
cfg.markersize = 3;
cfg.comment = 'no';
cfg.zlim = [-0.06 0.06];
cfg.xlim = [8 12];

set(groot,'defaultAxesFontName','Arial')
set(groot,'defaultTextFontName','Arial')

figure('Color','w','Position',[100 100 800 600]);
cfg.figure = gcf;
cfg.colormap = ft_colormap('balance',256);

ax(1) = subplot(2,2,1)
ft_topoplotER(cfg, GA_pre_with)
title({['\color[rgb]{' num2str(col_pre) '}Pre-Stimulus'], ...
       ['\color[rgb]{' num2str(col_self) '}Fatigue']}, ...
       'FontSize', 15, 'FontWeight', 'bold', 'Interpreter', 'tex')

ax(2) = subplot(2,2,2)
ft_topoplotER(cfg, GA_pre_without)
title({['\color[rgb]{' num2str(col_pre) '}Pre-Stimulus'], ...
       ['\color[rgb]{' num2str(col_no) '}No fatigue']}, ...
       'FontSize', 15, 'FontWeight', 'bold', 'Interpreter', 'tex')

ax(3) = subplot(2,2,3)
ft_topoplotER(cfg, GA_post_with)
title({['\color[rgb]{' num2str(col_post) '}Post-Stimulus'], ...
       ['\color[rgb]{' num2str(col_self) '}Fatigue']}, ...
       'FontSize', 15, 'FontWeight', 'bold', 'Interpreter', 'tex')

ax(4) = subplot(2,2,4)
ft_topoplotER(cfg, GA_post_without)
title({['\color[rgb]{' num2str(col_post) '}Post-Stimulus'], ...
       ['\color[rgb]{' num2str(col_no) '}No fatigue']}, ...
       'FontSize', 15, 'FontWeight', 'bold', 'Interpreter', 'tex')

ax(1).Position = [0.08 0.56 0.28 0.30];
ax(2).Position = [0.34 0.56 0.28 0.30];
ax(3).Position = [0.08 0.16 0.28 0.30];
ax(4).Position = [0.34 0.16 0.28 0.30];

letters = {'A','B','C','D'};

for k = 1:4
    axes(ax(k))
    text(-0.05,1.05,letters{k},...
        'Units','normalized',...
        'FontName','Arial',...
        'FontWeight','bold',...
        'FontSize',16,...
        'HorizontalAlignment','left');
end

% Add colorbar
h = colorbar('Position',[0.70 0.16 0.018 0.74]);

ylabel(h,'Power (\muV^2)',...
    'FontName','Arial',...
    'FontSize',16,...
    'FontWeight','bold');

set(h,...
    'FontName','Arial',...
    'FontSize',14);

% Export
set(gcf,'Units','pixels','Position',[100 100 2500 1800])

exportgraphics(gcf,...
'C:\Users\chris\Desktop\EPOC\EPOC-Kiel\03_figures\06_eeg_peaks\topo_alpha.pdf',...
'Resolution',600)
