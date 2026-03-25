clc; clear all; close all;

%% prepare everything

% path to folder with ICA weights
path                = "C:\Users\chris\Desktop\EPOC\01_data\01_individual_data_after_ICA";

% list of all .mat files
data = dir(fullfile(path, '*.mat'));

MAIN                = 'C:\Users\chris\Desktop\EPOC\';                       % paste path to EPOC folder
cd(MAIN);                                                                   % go into EPOC folder

epocPath            = genpath(MAIN);
fieldtripPath       = genpath('C:\Users\chris\Desktop\EPOC\04_software\fieldtrip-20250106');
epocPath            = strrep(epocPath, fieldtripPath, '');                  % add fieldtrip seperately from the rest 

addpath(epocPath);                                                          % add EPOC 
addpath('C:\Users\chris\Desktop\EPOC\04_software\fieldtrip-20250106');      % and fieldtrip 

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

    %% average pre
    cfg                 = [];
    cfg.latency         = [-1 0];
    data_pre            = ft_selectdata(cfg, data_erp);
    
    % average trials
    cfg                 = [];
    data_timelock       = ft_timelockanalysis(cfg, data_pre);

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

    % mean channels
    original.powspctrm          = mean(original.powspctrm, 1);         
    original.label              = {'avg_channels'};   
    fractal.powspctrm           = mean(fractal.powspctrm, 1);
    fractal.label               = {'avg_channels'};

    % single-trial powerspectrum
    cfg                         = [];
    cfg.output                  = 'pow';                                    % only powerspectrum
    cfg.method                  = 'mtmfft';
    cfg.tapsmofrq               = 4;
    cfg.pad                     = 4;
    cfg.foilim                  = [1 40]; 
    freq_single                 = ft_freqanalysis(cfg, data_pre);  
    % average over channels for trial
    pow                         = squeeze(mean(freq_single.powspctrm, 3));  % [nTrials x nFreqs]
    mean_pow                    = mean(pow, 1);                             % [1 x nFreqs]
    
    %% display the spectra 
    figure('WindowState', 'maximized');
    subplot(1,2,1); hold on;
    plot((original.freq), (original.powspctrm),'k');
    plot((fractal.freq), (fractal.powspctrm));
    xlabel('freq'); ylabel('power'); grid on;
    legend({'original','fractal'},'location','southwest');

    subplot(1,2,2); hold on;
    % single trial spectra
    for j = 1:size(pow, 1)
        plot(freq_single.freq, pow(j, :), 'Color', [0.85 0.85 0.85]);  
    end

    % average of all trials in red and thick 
    plot(freq_single.freq, mean_pow, 'r', 'LineWidth', 2);  

    xlabel('Frequency (Hz)');
    ylabel('Power (µV²)');
    title('Single-Trial Spectra + Mean');
    legend({'Single Trials', 'Mean'}, 'Location', 'northeast');
    grid on;

    saveas(gcf, ['C:\Users\chris\Desktop\EPOC\03_figures\03_fooof_pvt' filesep sub '_prestim.png']); 
    close;

    % average and save aperiodic data
    aperiodic_all           = fractal.fooofparams.aperiodic_params;         % [nChan x 3]
    r_squared_all           = fractal.fooofparams.r_squared;                % [nChan x 1]

    offset                  = mean(aperiodic_all(:,1));
    knee                    = mean(aperiodic_all(:,2));
    exponent                = mean(aperiodic_all(:,3));
    rsq                     = mean(r_squared_all);

    results_pre(i).id       = sub;
    results_pre(i).offset   = offset;
    results_pre(i).knee     = knee;
    results_pre(i).exponent = exponent;
    results_pre(i).rsq      = rsq;

    % add alpha
    corr_power              = original.powspctrm - fractal.powspctrm; 
    % Alpha-Peak 
    mask_alpha              = original.freq >= 8 & original.freq <= 12;
    [~, max_idx]            = max(corr_power(mask_alpha));
    alpha_peak_freqs        = original.freq(mask_alpha);
    alpha_peak              = alpha_peak_freqs(max_idx);
    alpha_peak_power        = corr_power(original.freq == alpha_peak);
    mean_alpha_power        = mean(corr_power(mask_alpha));

    % show alpha
    figure; hold on;
    plot(original.freq, corr_power, 'b');
    plot(alpha_peak, corr_power(original.freq==alpha_peak), 'ro', 'MarkerSize',8,'LineWidth',2);
    xlabel('Frequency (Hz)'); ylabel('Power (µV², corr)');
    title(['Alpha Peak: ' num2str(alpha_peak, '%.2f') ' Hz']);
    xline(8, '--k', 'LineWidth', 0.8);
    xline(12, '--k', 'LineWidth', 0.8);
    grid on;
    saveas(gcf, ['C:\Users\chris\Desktop\EPOC\03_figures\04_alpha_peak_pvt' filesep sub '_prestim.png']);
    close;

    % add theta
    corr_power              = original.powspctrm - fractal.powspctrm; 
    % Theta-Peak 
    mask_theta              = original.freq >= 4 & original.freq <= 7;
    [~, max_idx]            = max(corr_power(mask_theta));
    theta_peak_freqs        = original.freq(mask_theta);
    theta_peak              = theta_peak_freqs(max_idx);
    theta_peak_power        = corr_power(original.freq == theta_peak);
    mean_theta_power        = mean(corr_power(mask_theta));

    % show alpha
    figure; hold on;
    plot(original.freq, corr_power, 'b');
    plot(theta_peak, corr_power(original.freq==theta_peak), 'ro', 'MarkerSize',8,'LineWidth',2);
    xlabel('Frequency (Hz)'); ylabel('Power (µV², corr)');
    title(['Theta Peak: ' num2str(theta_peak, '%.2f') ' Hz']);
    xline(4, '--k', 'LineWidth', 0.8);
    xline(7, '--k', 'LineWidth', 0.8);
    grid on;
    saveas(gcf, ['C:\Users\chris\Desktop\EPOC\03_figures\05_theta_peak_pvt' filesep sub '_prestim.png']);
    close;

    % speichern
    results_pre(i).alpha_peak_freq      = alpha_peak;
    results_pre(i).alpha_peak_power     = alpha_peak_power;
    results_pre(i).mean_alpha_power     = mean_alpha_power;

    results_pre(i).theta_peak_freq      = theta_peak;
    results_pre(i).theta_peak_power     = theta_peak_power;
    results_pre(i).mean_theta_power     = mean_theta_power;

    %% average post

    cfg             = [];
    cfg.latency     = [0 1];
    data_post       = ft_selectdata(cfg, data_erp);
    
    % average trials
    cfg             = [];
    data_timelock   = ft_timelockanalysis(cfg, data_post);

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

    % mean channels
    original.powspctrm          = mean(original.powspctrm, 1);         
    original.label              = {'avg_channels'};   
    fractal.powspctrm           = mean(fractal.powspctrm, 1);
    fractal.label               = {'avg_channels'};

    % single-trial fooof
    cfg                         = [];
    cfg.output                  = 'pow';                                    % only powerspectrum
    cfg.method                  = 'mtmfft';
    cfg.tapsmofrq               = 4;
    cfg.pad                     = 4;
    cfg.foilim                  = [1 40]; 
    freq_single                 = ft_freqanalysis(cfg, data_post);  
    % average over channels for trial
    pow                         = squeeze(mean(freq_single.powspctrm, 3));  % [nTrials x nFreqs]
    mean_pow                    = mean(pow, 1);                             % [1 x nFreqs]

    %% display the spectra 
    figure('WindowState', 'maximized');
    subplot(1,2,1); hold on;
    plot((original.freq), (original.powspctrm),'k');
    plot((fractal.freq), (fractal.powspctrm));
    xlabel('freq'); ylabel('power'); grid on;
    legend({'original','fractal'},'location','southwest');

    subplot(1,2,2); hold on;
    % single trial spectra
    for j = 1:size(pow, 1)
        plot(freq_single.freq, pow(j, :), 'Color', [0.85 0.85 0.85]);  
    end

    % average of all trials in red and thick 
    plot(freq_single.freq, mean_pow, 'r', 'LineWidth', 2);  

    xlabel('Frequency (Hz)');
    ylabel('Power (µV²)');
    title('Single-Trial Spectra + Mean');
    legend({'Single Trials', 'Mean'}, 'Location', 'northeast');
    grid on;

    saveas(gcf, ['C:\Users\chris\Desktop\EPOC\03_figures\03_fooof_pvt' filesep sub '_poststim.png']); 
    close;

    % average and save aperiodic data
    aperiodic_all               = fractal.fooofparams.aperiodic_params;     % [nChan x 3]
    r_squared_all               = fractal.fooofparams.r_squared;            % [nChan x 1]

    offset                      = mean(aperiodic_all(:,1));
    knee                        = mean(aperiodic_all(:,2));
    exponent                    = mean(aperiodic_all(:,3));
    rsq                         = mean(r_squared_all);

    results_post(i).id          = sub;
    results_post(i).offset      = offset;
    results_post(i).knee        = knee;
    results_post(i).exponent    = exponent;
    results_post(i).rsq         = rsq;

    % add alpha
    corr_power                  = original.powspctrm - fractal.powspctrm; 
    % Alpha-Peak 
    mask_alpha                  = original.freq >= 8 & original.freq <= 12;
    [~, max_idx]                = max(corr_power(mask_alpha));
    alpha_peak_freqs            = original.freq(mask_alpha);
    alpha_peak                  = alpha_peak_freqs(max_idx);
    alpha_peak_power            = corr_power(original.freq == alpha_peak);
    mean_alpha_power            = mean(corr_power(mask_alpha));

    % show alpha
    figure; hold on;
    plot(original.freq, corr_power, 'b');
    plot(alpha_peak, corr_power(original.freq==alpha_peak), 'ro', 'MarkerSize',8,'LineWidth',2);
    xlabel('Frequency (Hz)'); ylabel('Power (µV², corr)');
    title(['Alpha Peak: ' num2str(alpha_peak, '%.2f') ' Hz']);
    xline(8, '--k', 'LineWidth', 0.8);
    xline(12, '--k', 'LineWidth', 0.8);
    grid on;
    saveas(gcf, ['C:\Users\chris\Desktop\EPOC\03_figures\04_alpha_peak_pvt' filesep sub '_poststim.png']);
    close;

    % add theta
    corr_power                  = original.powspctrm - fractal.powspctrm; 
    % theta-Peak 
    mask_theta                  = original.freq >= 4 & original.freq <= 7;
    [~, max_idx]                = max(corr_power(mask_theta));
    theta_peak_freqs            = original.freq(mask_theta);
    theta_peak                  = theta_peak_freqs(max_idx);
    theta_peak_power            = corr_power(original.freq == theta_peak);
    mean_theta_power            = mean(corr_power(mask_theta));

    % show theta
    figure; hold on;
    plot(original.freq, corr_power, 'b');
    plot(theta_peak, corr_power(original.freq==theta_peak), 'ro', 'MarkerSize',8,'LineWidth',2);
    xlabel('Frequency (Hz)'); ylabel('Power (µV², corr)');
    title(['theta Peak: ' num2str(theta_peak, '%.2f') ' Hz']);
    xline(4, '--k', 'LineWidth', 0.8);
    xline(7, '--k', 'LineWidth', 0.8);
    grid on;
    saveas(gcf, ['C:\Users\chris\Desktop\EPOC\03_figures\05_theta_peak_pvt' filesep sub '_poststim.png']);
    close;

    % speichern
    results_post(i).alpha_peak_freq     = alpha_peak;
    results_post(i).alpha_peak_power    = alpha_peak_power;
    results_post(i).mean_alpha_power    = mean_alpha_power;

    results_post(i).theta_peak_freq     = theta_peak;
    results_post(i).theta_peak_power    = theta_peak_power;
    results_post(i).mean_theta_power    = mean_theta_power;

end

%%
save([PATHOUT '02_fooof_pvt' filesep 'fooofed.mat'], "results_pre", "results_post");
save([PATHOUT 'bad_trials_fooof_s02.mat'], "all_bad_trials");