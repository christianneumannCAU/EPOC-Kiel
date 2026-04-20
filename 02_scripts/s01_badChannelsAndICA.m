% clear command window, workspace and figures
clc; clear all; close all;

%% Step 0: add + prepare everything

MAIN                = 'C:\Users\chris\Desktop\EPOC\EPOC-Kiel\';                       % paste path to EPOC folder
cd(MAIN);                                                                   % go into EPOC folder

epocPath            = genpath(MAIN);
fieldtripPath       = genpath('C:\Users\chris\Desktop\EPOC\EPOC-Kiel\04_software\fieldtrip-20250106');
epocPath            = strrep(epocPath, fieldtripPath, '');                  % add fieldtrip seperately from the rest 

addpath(epocPath);                                                          % add EPOC 
addpath('C:\Users\chris\Desktop\EPOC\EPOC-Kiel\04_software\fieldtrip-20250106');      % and fieldtrip 

% Set envir
PATHIN              = [MAIN '01_data\00_bids\'];                            % Path to BIDS data
PATHOUT             = [MAIN '01_data\'];                                    % Path where to save the data

% prepare layout
cfg = [];
load('BC-128-pass-lay.mat');
cfg.elec            = elec;
layout              = ft_prepare_layout(cfg);

% prepare patients 
patient             = dir(PATHIN);
patient([1,2],:)    = [];

isVP                = startsWith({patient.name}, "sub");
patient             = patient(isVP);

% prepare ICA table
all_IC_tables       = [];

% rmove path to EEG lab
rmpath(genpath('C:\Users\chris\Desktop\EPOC\EPOC-Kiel\04_software\eeglab_current\eeglab2024.2')); 
rehash;
% initiate fieldtrip
ft_defaults;
which runica  

%% Loop over every patient

for p               = 1 : numel(patient)
    %% Step 1: read data
    
    cd([PATHIN patient(p).name filesep 'eeg']);                             % go to current patient
    % search for pvt data
    indat           = dir(fullfile(PATHIN, patient(p).name, 'eeg', '*pvt_eeg.vhdr'));      

    % check, if indat empty
    if isempty(indat)
        disp(['Keine passende Datei gefunden in: ' fullfile(PATHIN, patient(p).name, 'eeg')]);
        continue;  % skips this Iteration of the Loop
    end

    % check, if data for this VP already exists
    exists          = dir('C:\Users\chris\Desktop\EPOC\EPOC-Kiel\01_data\01_individual_data_after_ICA_pvt');
    ex              = cellfun(@(x) extractBefore(x, '.'), {exists.name}, 'UniformOutput', false);
    if ismember(extractBefore(indat.name, 11), ex) 
        continue;
    end

    % define trial-definition and save for ICA 
    cfg                 = [];
    cfg.dataset         = [PATHIN patient(p).name filesep 'eeg' filesep indat.name];
    cfg.trialfun        = 'ft_trialfun_general';
    cfg.trialdef.length = 1;                                                % cut data into 1s segments
    cfg                 = ft_definetrial(cfg); 
    trl                 = cfg.trl;                                          % save the trial-definition

    % read continous raw data
    cfg = [];
    cfg.dataset         = [PATHIN patient(p).name filesep 'eeg' filesep indat.name];
    cfg.continous       = 'yes'; 
    
    data_r              = ft_preprocessing(cfg);                            % save raw data!
    
    %% Change channel names

    for c           = 1:length(data_r.label)
        data_r.label{c} = layout.label{c};
    end
    
    %% Step 2: preprocessing for ica 

    cfg = [];
    cfg.demean      = 'yes';                                                % remove DC offset
    cfg.hpfilter    = 'yes';
    cfg.hpfreq      = 1;                                                    % high-pass filter, cutting everything under 1 Hz
    cfg.hpfilttype  = 'firws';
    cfg.lpfilter    = 'yes';
    cfg.lpfreq      = 40;                                                   % low-pass filter, cutting everything over 40 Hz 
    cfg.lpfilttype  = 'firws';
    cfg.pad         = 'nextpow2';

    data_pica       = ft_preprocessing(cfg, data_r);

    %% Step 3: reject bad channels

    % calculate variance for every channel
    channel_variance    = var(data_pica.trial{1}, 0, 2);                     
    
    median_var          = median(channel_variance);
    mad_var             = mad(channel_variance, 1);                         % Median Absolute Deviation
    adaptive_threshold  = median_var + 1 * mad_var;                         % define threshold 

    cfg                 = [];
    cfg.metric          = 'var';
    if adaptive_threshold > 1500
        cfg.threshold   = 1500;                                             % add absolute criterium
    else
        cfg.threshold   = adaptive_threshold;
    end
    bad_channels        = ft_badchannel(cfg, data_pica);

    % only keep good channels
    cfg                 = [];
    cfg.channel         = setdiff(data_pica.label, bad_channels.badchannel);       
    data_c              = ft_selectdata(cfg, data_pica);

    % save how many channels were removed per participant
    bad_channels_all{p} = bad_channels.badchannel;

    %% Step 4: ICA Trial definition 

    % Cut the data according to the trial definition from the beginning
    cfg             = [];
    cfg.trl         = trl;                                                  
    data_tica       = ft_redefinetrial(cfg,data_c);

    %% Step 5: remove bad segments for ICA 

    % calculate mean for every segment across all channels and samples
    trial_means     = cellfun(@(x) mean(x(:)), data_tica.trial);

    % calculate z-score for every segment
    z_scores        = (trial_means - mean(trial_means)) / std(trial_means);

    % define threshold 
    threshold       = 2;

    % identify bad segments
    bad_trials      = find(abs(z_scores) > threshold);

    % delete bad segments
    cfg             = [];
    cfg.trials      = setdiff(1:length(data_tica.trial), bad_trials);
    data_tica       = ft_selectdata(cfg, data_tica);

    % plot all z-scores
    scatter(1:length(z_scores), z_scores, 'filled');
    xlabel('Trial');
    ylabel('Z-Score');
    title('Z-Score per Trial');
    yline(2, 'r--', 'Threshold +2');
    yline(-2, 'r--', 'Threshold -2')
    saveas(gcf, ['C:\Users\chris\Desktop\EPOC\EPOC-Kiel\03_figures\01_bad_trials_z_pvt' filesep patient(p).name '_trials_z.png']); 
    close;

    % save how many segments were removed per participant
    bad_segments_all{p} = bad_trials;

    %% Step 6: automatic ICA
    
    % add path to EEGLAB
    addpath(genpath('C:\Users\chris\Desktop\EPOC\EPOC-Kiel\04_software\eeglab_current\eeglab2024.2'));
    eeglab; close;                                                          
    
    EEG                     = fieldtrip2eeglab(data_tica, data_tica.trial); % change fieldtrip data to eeglab data
    EEG                     = eeg_checkset(EEG);

    try
        EEG                 = pop_runica(EEG, 'icatype', 'runica');         % perform ICA
    catch
        disp(['Skipping participant error at runica: ', patient(p).name]);
        continue;  
    end

    EEG.srate               = round(EEG.srate);                             % round Sampling-Rate 

    % add the electrode position for the remaining channels
    [~, idx]                = ismember({EEG.urchanlocs.labels}, elec.label);
    elec_filtered.chanpos   = elec.chanpos(idx, :);
    elec_filtered.label     = elec.label(idx);
    EEG.chanlocs            = struct('labels', elec_filtered.label, ...
                          'X', num2cell(elec_filtered.chanpos(:,1)), ...
                          'Y', num2cell(elec_filtered.chanpos(:,2)), ...
                          'Z', num2cell(elec_filtered.chanpos(:,3)));

    EEG                     = iclabel(EEG);                                 % label components

    % reject components (see function help message)
    EEG                     = pop_icflag(EEG, [0 0;0.3 1; 0.8 1; 0.4 1; 0.3 1; 0.3 1; 0.35 1]); 
    rejected_comps          = find(EEG.reject.gcompreject > 0);       
    % what are the good components (for plot)
    good_comps              = setdiff(1:size(EEG.icawinv,2), rejected_comps);
    % save rejeced components 
    rejected_comps_all{p}   = rejected_comps;                                 

    
    % create summarizing table
    % table for one participant
    IC_numbers          = (1:size(EEG.etc.ic_classification.ICLabel.classifications, 1))';  
    IC_labels           = EEG.etc.ic_classification.ICLabel.classes;  
    IC_probs            = EEG.etc.ic_classification.ICLabel.classifications;  
    IC_table            = array2table(IC_probs, 'VariableNames', IC_labels);  
    IC_table.IC_Number  = IC_numbers; 

    % table for all participants
    ID_column           = repmat({patient(p).name}, size(IC_table,1), 1);
    [max_prob, max_idx] = max(table2array(IC_table(:,1:7)), [], 2);
    max_label           = IC_table.Properties.VariableNames(max_idx)';

    current_table       = table(ID_column, IC_table.IC_Number, max_label, max_prob, ...
        'VariableNames', {'ID', 'IC_Number', 'Max_Label', 'Max_Probability'});
    
    if isempty(all_IC_tables)
        all_IC_tables   = current_table;
    else
        all_IC_tables   = [all_IC_tables; current_table];
    end

    % create comp structure for fieldtrip
    comp            = [];
    comp.label      = strcat('IC', string(1:size(EEG.icaweights,1)))';      % IC-names
    comp.topolabel  = {EEG.chanlocs.labels};                                % original EEG channels
    comp.unmixing   = EEG.icaweights * EEG.icasphere;                       % unmixing matrix
    comp.topo       = EEG.icawinv;                                          % mixed matric (inverse ICA)
    comp.time       = EEG.times / 1000;                                     % time in seconds
    
    if isempty(EEG.icaact)
        EEG.icaact  = (EEG.icaweights * EEG.icasphere) * EEG.data(EEG.icachansind, :);
    end
    
    % if epoched, reshape to trials
    if EEG.trials > 1
        EEG.icaact  = reshape(EEG.icaact, size(EEG.icaweights,1), EEG.pnts, EEG.trials);
        comp.trial  = arrayfun(@(x) EEG.icaact(:,:,x), 1:EEG.trials, 'UniformOutput', false);
    else
        comp.trial  = {EEG.icaact};                                         % if continous, save as cell 
    end
    
    %% Compare before & after
    
    % figure;plot(data_ci.time{1},data_ci.trial{5}(1,:),'r');hold;plot(data_tica.time{1},data_tica.trial{5}(1,:),'b')

    %% plot components
    figure('Name', sprintf('%s - Gute Komponenten', patient(p).name), 'NumberTitle', 'off');
    set(gcf, 'Position', get(0, 'ScreenSize'));
    cfg                 = [];
    cfg.component       = good_comps;

    load('C:\Users\chris\Desktop\EPOC\BC-128-pass-lay.mat');
    cfg.elec            = elec;
    layout              = ft_prepare_layout(cfg);
    cfg.layout          = layout;

    % good components
    cfg.comment         = 'no'
    ft_topoplotIC(cfg, comp);
    saveas(gcf, ['C:\Users\chris\Desktop\EPOC\EPOC-Kiel\03_figures\02_ica_plots_pvt' filesep patient(p).name '_goodcomps.png']); 
    close;

    % bad components
    figure('Name', sprintf('%s - Entfernte Komponenten', patient(p).name), 'NumberTitle', 'off');
    set(gcf, 'Position', get(0, 'ScreenSize'));
    cfg.component       = rejected_comps;
    ft_topoplotIC(cfg, comp);
    saveas(gcf, ['C:\Users\chris\Desktop\EPOC\03_figures\02_ica_plots_pvt' filesep patient(p).name '_badcomps.png']); 
    close;

    %% save data
    
    save([PATHOUT '01_individual_data_after_ICA' filesep patient(p).name '.mat'], "comp", "data_r", "bad_channels", "bad_channels_all", "IC_table", "rejected_comps");

    %% clean up before next loop
    clearvars -except MAIN PATHIN PATHOUT layout p patient all_IC_tables bad_channels_all bad_segments_all rejected_comps_all elec  

end

save([PATHOUT 'noisy_objects_after_s01.mat'], "bad_channels_all", "bad_segments_all", "rejected_comps_all", "all_IC_tables");