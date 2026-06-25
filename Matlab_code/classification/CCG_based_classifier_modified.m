function CCG_based_classifier_modified(custom_settings)
% Generate CCG review sessions for curated excitatory/inhibitory labels.

if nargin < 1 || isempty(custom_settings)
    custom_settings = struct();
end

config = load_classification_config(get_override_value(custom_settings, 'configPath', ''));
add_repo_matlab_code_path();
add_configured_dependency_paths(config, custom_settings);

make_plots = logical(get_override_value(custom_settings, 'makePlots', true));

sessionInfoPath = char(string(get_override_value(custom_settings, 'sessionInfoPath', '')));
if isempty(sessionInfoPath)
    sessionInfoPath = resolve_session_info_path();
end
loadedSessionInfo = load(sessionInfoPath, 'sessInfo');
if ~isfield(loadedSessionInfo, 'sessInfo')
    error('sessInfo was not found in %s.', sessionInfoPath);
end
sessInfo = loadedSessionInfo.sessInfo;

allCellsPath = char(string(get_override_value(custom_settings, 'allCellsPath', '')));
if isempty(allCellsPath)
    allCellsPath = resolve_all_cells_path();
end
loadedAllCells = load(allCellsPath, 'All_Cells_combined');
if ~isfield(loadedAllCells, 'All_Cells_combined')
    error('All_Cells_combined was not found in %s.', allCellsPath);
end
All_Cells_combined = loadedAllCells.All_Cells_combined;

review_root = char(string(get_override_value(custom_settings, 'reviewRoot', '')));
if isempty(review_root)
    review_root = default_review_root();
end
if ~exist(review_root, 'dir')
    mkdir(review_root)
end

review_generation_mode = char(string(get_override_value(custom_settings, 'reviewGenerationMode', '')));
if isempty(review_generation_mode)
    review_generation_mode = resolve_review_generation_mode('only_missing_review_files');
end
if isempty(review_generation_mode)
    fprintf('CCG classifier run cancelled before processing.\n');
    return
end

group1 = get_override_value(custom_settings, 'sessionIndices', default_ccg_session_indices());

ttFiles_matrix = [];
processed_session_count = 0;
skipped_existing_review_count = 0;

for i_iter = 1:numel(group1)
    iii = group1(i_iter);
    review_file = fullfile(review_root, sprintf('Session%03d-Animal%s-Day%s-CCGReview.mat', ...
        iii, num2str(sessInfo(iii).animal), num2str(sessInfo(iii).day)));

    if should_skip_review_session(review_file, review_generation_mode)
        fprintf('Skipping existing review file for session %d: %s \n', iii, review_file);
        skipped_existing_review_count = skipped_existing_review_count + 1;
        continue
    end

    fprintf('Processing: %s \n', sessInfo(iii).mainDir);
    clear spikeData spikeTimes

    fprintf("Processing i: %s \n", num2str(iii));
    %% Read TT List
    clear tt_files
    cell_no = 0;
    fid=fopen(fullfile(sessInfo(iii).mainDir, sessInfo(iii).tList));
    while 1
        cell_no = cell_no +1;
        tline = fgetl(fid);
        if ~ischar(tline), break, end
        tt_files{cell_no} = tline;
    end
    fclose(fid);
    numCells = 0;
    numCells = numel(tt_files);
    ttFiles_matrix = [ttFiles_matrix, [tt_files]];

    
    [number_OFs , folder_numbers] = count_of_folders(fullfile(sessInfo(iii).mainDir));
    behavior_spike_times = cell(number_OFs,1);
    for sessions = 1:number_OFs
        clear spikeTimes spikeData
        fprintf(' - - Reading Spikes %s :  ', sessInfo(iii).sessDirs{sessions}); 
        % Read spikes
        spikeData = readSpikeDataOnly(fullfile(sessInfo(iii).mainDir, sessInfo(iii).sessDirs{sessions}) , tt_files);
        spikeTimes = fixSpikes(spikeData);
        behavior_spike_times{sessions} = spikeTimes;
        if sessions == 1 % initialize
            allOF_tSp = cell(size(spikeTimes,1),1);
        end
        
        session_time(sessions) = range(vertcat(spikeTimes{:})); % calculate once the total recording time of each OF
        
        for ccc = 1:size(spikeTimes,1)
            clear tSp
            allOF_tSp{ccc} = [allOF_tSp{ccc}; spikeTimes{ccc}];
            
            % Need to sum here duration of recordings, to normalize to Hz
            
            % [intrins_frequ_sess(ccc,sessions),thetaMod_score_sess(ccc,sessions), ac_timing, ac_nSp]  = ACG_ThetaMod(tSp,'none', 0);
        end
         clear spikeTimes spikeData
    end
 
    for sleeps = 1:2
        clear spikeTimes spikeData
        fprintf(' - - Reading Spikes %s: ', sessInfo(iii).sleepDirs{sleeps}); 
        spikeData = readSpikeDataOnly(fullfile(sessInfo(iii).mainDir, sessInfo(iii).sleepDirs{sleeps}) , tt_files);
        spikeTimes = fixSpikes(spikeData);
        session_time(sleeps+3) = range(vertcat(spikeTimes{:}));
        for ccc = 1:size(spikeTimes,1)
            clear tSp
            allOF_tSp{ccc} = [allOF_tSp{ccc}; spikeTimes{ccc}];
        end
    end

    clear spikeTimes spikeData
    
    spikeTimes = allOF_tSp;

clear spikes

for i = 1: size(spikeTimes,1)
num_spikes =   size(spikeTimes{i,1},1);
tt_cell = tt_files{i};
tok = regexp(tt_cell, '^TT(\d+)', 'tokens', 'once');
if isempty(tok)
    error('Could not parse tetrode ID from tt file name: %s', tt_cell);
end
tetrode_id = str2double(tok{1});

spikes.shankID(i,1) = tetrode_id;
%spikes.shankID(i,1) = i;
spikes.cluID(i,1) = i;

end
spikes.times  = allOF_tSp;
cell_classification = zeros(size(allOF_tSp,1), 1);
%mono_res = ce_MonoSynConvClick(spikes,'includeInhibitoryConnections',true); % detects the monosynaptic connections

%try

% Keep original spikes struct
spikes_full = spikes;
nUnits_full = numel(spikes_full.times);

% Units that actually have spikes in this behavioral block
valid_units = ~cellfun(@isempty, spikes_full.times);
[total_pairs_evaluated, total_pairs_crossTT, total_pairs_sameTT] = ...
    compute_total_pair_counts(spikes_full.shankID, valid_units);

sig_window_exc_ms = [1.0 2.8];
sig_window_inh_ms = [1.0 4.0];
plot_bin_full_s = 0.0010;      % 1.0 ms for full plot
plot_bin_zoom_exc_s = 0.0002;  % 0.2 ms detector-matched excitatory zoom
plot_bin_zoom_inh_s = 0.0002;  % 0.2 ms detector-matched inhibitory zoom
plot_bin_raw_s = 0.0001;       % 0.1 ms raw display-only plot
zoom_window_exc_ms = [-8 8];
zoom_window_inh_ms = [-8 8];
save_inhibitory_residual_plots = true;

detection_config = struct( ...
    'bin_s', 0.0002, ...
    'duration_s', 0.120, ...
    'sig_window_exc_s', 0.0028, ...
    'sig_window_inh_s', 0.0040, ...
    'alpha', 0.001, ...
    'conv_w_s', 0.010, ...
    'reference_window_s', [-0.004 -0.001], ...
    'causal_window_start_s', 0.001, ...
    'same_shank_mask_s', [-0.001 0.001], ...
    'rebound_window_s', [0.001 0.004], ...
    'include_inhibitory_connections', true, ...
    'sr', 32000 ...
    );
settings_signature = make_settings_signature(detection_config);

if sum(valid_units) < 2
    warning('Fewer than 2 non-empty units in this session. Skipping ce_MonoSynConvClick.');
    mono_res_excitatory = [];
    mono_res_inhibitory = [];
    
    cell_classification = zeros(nUnits_full,1);
    review_session = struct();
    review_session.version = 2;
    review_session.session_index = iii;
    review_session.session_label = sprintf('Animal %s Day %s Session %s', ...
        num2str(sessInfo(iii).animal), num2str(sessInfo(iii).day), num2str(iii));
    review_session.mainDir = sessInfo(iii).mainDir;
    review_session.animal = sessInfo(iii).animal;
    review_session.day = sessInfo(iii).day;
    review_session.behavior_session_dirs = sessInfo(iii).sessDirs(1:number_OFs);
    review_session.n_units = nUnits_full;
    review_session.tt_files = tt_files;
    review_session.firing_rates = All_Cells_combined(iii).classific_firingRate;
    review_session.auto_cell_classification = cell_classification;
    review_session.behavior_spike_times = behavior_spike_times;
    review_session.exc_pairs_all = zeros(0,2);
    review_session.inh_pairs_all = zeros(0,2);
    review_session.rebound_override_pairs = zeros(0,2);
    review_session.review_pairs = initialize_review_pairs(zeros(0,2), zeros(0,2), zeros(0,2), spikes_full.shankID);
    review_session.spike_times = spikes_full.times;
    review_session.mono_res_excitatory = mono_res_excitatory;
    review_session.mono_res_inhibitory = mono_res_inhibitory;
    review_session.plot_config = struct( ...
        'sig_window_exc_ms', sig_window_exc_ms, ...
        'sig_window_inh_ms', sig_window_inh_ms, ...
        'plot_bin_full_s', plot_bin_full_s, ...
        'plot_bin_zoom_exc_s', plot_bin_zoom_exc_s, ...
        'plot_bin_zoom_inh_s', plot_bin_zoom_inh_s, ...
        'plot_bin_raw_s', plot_bin_raw_s, ...
        'zoom_window_exc_ms', zoom_window_exc_ms, ...
        'zoom_window_inh_ms', zoom_window_inh_ms, ...
        'full_window_ms', [-60 60] ...
        );
    review_session.detection_config = detection_config;
    review_session.settings_signature = settings_signature;
    review_session.auto_pair_counts = struct( ...
        'total_connections_tested', total_pairs_evaluated, ...
        'total_pairs_evaluated', total_pairs_evaluated, ...
        'total_pairs_crossTT', total_pairs_crossTT, ...
        'total_pairs_sameTT', total_pairs_sameTT, ...
        'identified_candidate_pairs', 0, ...
        'detected_exc_all', 0, ...
        'detected_inh_all', 0, ...
        'detected_exc_crossTT', 0, ...
        'detected_exc_sameTT', 0, ...
        'detected_inh_crossTT', 0, ...
        'detected_inh_sameTT', 0 ...
        );

    review_session = merge_review_session_state(review_session, review_file);
    review_session.classifier_run_mode = review_generation_mode;
    save(review_file, 'review_session', '-v7.3');
    fprintf('Saved skipped review session to %s\n', review_file);
    processed_session_count = processed_session_count + 1;
    continue
end

% Reduced spikes struct for ce_MonoSynConvClick
spikes = struct();
spikes.times   = spikes_full.times(valid_units);
spikes.shankID = spikes_full.shankID(valid_units,:);
spikes.cluID   = spikes_full.cluID(valid_units,:);

sr = detection_config.sr;
mono_res_excitatory = ce_MonoSynConvClick(spikes, ...
    'binsize', detection_config.bin_s, ... % 0.2 ms detection bin
    'duration', detection_config.duration_s, ...  % 120 ms CCG window
    'sigWindow', detection_config.sig_window_exc_s, ... % Causal window upper edge: 1.0-2.8 ms in this implementation
    'alpha', detection_config.alpha, ...  % Conservative significance threshold
    'conv_w', detection_config.conv_w_s, ...  % 10 ms convolution window
    'includeInhibitoryConnections', true, ...
    'sr', sr);

mono_res_inhibitory = ce_MonoSynConvClick(spikes, ...
    'binsize', detection_config.bin_s, ... % 0.2 ms detection bin
    'duration', detection_config.duration_s, ...  % 120 ms CCG window
    'sigWindow', detection_config.sig_window_inh_s, ... % Causal window upper edge: 1.0-4.0 ms in this implementation
    'alpha', detection_config.alpha, ...  % Conservative significance threshold
    'conv_w', detection_config.conv_w_s, ...  % 10 ms convolution window
    'includeInhibitoryConnections', true, ...
    'sr', sr);
mono_res_excitatory = expand_mono_res_to_full_indexing(mono_res_excitatory, valid_units, nUnits_full);
mono_res_inhibitory = expand_mono_res_to_full_indexing(mono_res_inhibitory, valid_units, nUnits_full);

if ~isempty(mono_res_excitatory.sig_con_excitatory)
    assert(size(mono_res_excitatory.sig_con_excitatory,2) == 2, ...
        'sig_con_excitatory is not N x 2 after expansion');
end
if ~isempty(mono_res_inhibitory.sig_con_inhibitory)
    assert(size(mono_res_inhibitory.sig_con_inhibitory,2) == 2, ...
        'sig_con_inhibitory is not N x 2 after expansion');
end

ccg_ts_ms = mono_res_excitatory.ccgTs * 1000;

% Post-filter excitatory pairs that show an earlier inhibitory trough
% followed by a later positive rebound in the +1 to +4 ms causal window.
exc_pairs_raw = mono_res_excitatory.sig_con_excitatory;
inh_pairs_raw = mono_res_inhibitory.sig_con_inhibitory;
rebound_override_pairs = detect_inhibitory_rebound_pairs(mono_res_inhibitory, exc_pairs_raw);

exc_pairs_all = remove_pair_rows(exc_pairs_raw, rebound_override_pairs);
inh_pairs_all = unique_pair_rows([inh_pairs_raw; rebound_override_pairs]);

% Split post-filtered pairs by tetrode relationship
[exc_pairs_crossTT, exc_pairs_sameTT] = split_pairs_by_tetrode(exc_pairs_all, spikes_full.shankID);
[inh_pairs_crossTT, inh_pairs_sameTT] = split_pairs_by_tetrode(inh_pairs_all, spikes_full.shankID);

%% New cell classification
% Final high-confidence training labels:
% - excitatory: at least one excitatory connection and no inhibitory connections
% - inhibitory: at least one inhibitory connection and no excitatory connections
% - inhibitory evidence includes same- and cross-TT pairs

cell_classification = zeros(numel(mono_res_excitatory.n), 1);

exc_cells = [];
inh_cells_all = [];

if ~isempty(exc_pairs_all)
    exc_cells = exc_pairs_all(:,1);
end

if ~isempty(inh_pairs_all)
    inh_cells_all = inh_pairs_all(:,1);
end

all_cells = unique([exc_cells; inh_cells_all]);

for i = 1:length(all_cells)
    cell_id = all_cells(i);

    exc_count = sum(exc_cells == cell_id);
    inh_count = sum(inh_cells_all == cell_id);

    if exc_count > 0 && inh_count == 0
        cell_classification(cell_id) = 1;   % Excitatory
    elseif inh_count > 0 && exc_count == 0
        cell_classification(cell_id) = 2;   % Inhibitory
    else
        cell_classification(cell_id) = 3;   % Ambiguous / mixed
    end
end

accepted_exc_pairs_crossTT = zeros(0,2);
accepted_exc_pairs_sameTT = zeros(0,2);
accepted_inh_pairs_crossTT = zeros(0,2);
accepted_inh_pairs_sameTT = zeros(0,2);

if ~isempty(exc_pairs_crossTT)
    exc_keep_crossTT = cell_classification(exc_pairs_crossTT(:,1)) == 1;
    accepted_exc_pairs_crossTT = exc_pairs_crossTT(exc_keep_crossTT, :);
end
if ~isempty(exc_pairs_sameTT)
    exc_keep_sameTT = cell_classification(exc_pairs_sameTT(:,1)) == 1;
    accepted_exc_pairs_sameTT = exc_pairs_sameTT(exc_keep_sameTT, :);
end
if ~isempty(inh_pairs_crossTT)
    inh_keep_crossTT = cell_classification(inh_pairs_crossTT(:,1)) == 2;
    accepted_inh_pairs_crossTT = inh_pairs_crossTT(inh_keep_crossTT, :);
end
if ~isempty(inh_pairs_sameTT)
    inh_keep_sameTT = cell_classification(inh_pairs_sameTT(:,1)) == 2;
    accepted_inh_pairs_sameTT = inh_pairs_sameTT(inh_keep_sameTT, :);
end

fprintf('High-confidence excitatory cells: %s\n', num2str(sum(cell_classification==1)));
disp(find(cell_classification==1))

fprintf('High-confidence inhibitory cells: %s\n', num2str(sum(cell_classification==2)));
disp(find(cell_classification==2))

fprintf('High-confidence unclear cells: %s\n', num2str(sum(cell_classification==3)));
disp(find(cell_classification==3))

fprintf('Reassigned excitatory-to-inhibitory rebound pairs: %d\n', size(rebound_override_pairs,1));
fprintf('Detected excitatory pairs after rebound filter: total %d | cross-TT %d | same-TT %d\n', ...
    size(exc_pairs_all,1), size(exc_pairs_crossTT,1), size(exc_pairs_sameTT,1));
fprintf('Detected inhibitory pairs after rebound filter: total %d | cross-TT %d | same-TT %d\n', ...
    size(inh_pairs_all,1), size(inh_pairs_crossTT,1), size(inh_pairs_sameTT,1));
fprintf('Accepted same-TT inhibitory pairs for training labels: %d\n', size(accepted_inh_pairs_sameTT,1));
fprintf('Accepted cross-TT inhibitory pairs for training labels: %d\n', size(accepted_inh_pairs_crossTT,1));

%% Visualize connectivity matrix
connMatrix = zeros(size(mono_res_excitatory.n,1));

for k = 1:size(accepted_exc_pairs_crossTT, 1)
    pre = accepted_exc_pairs_crossTT(k, 1);
    post = accepted_exc_pairs_crossTT(k, 2);
    connMatrix(post, pre) = 1;
end

for k = 1:size(accepted_inh_pairs_crossTT, 1)
    pre = accepted_inh_pairs_crossTT(k, 1);
    post = accepted_inh_pairs_crossTT(k, 2);
    connMatrix(post, pre) = -1;
end

connMatrixPlot = connMatrix + 2;

cmap = [
    1 0 0;  % Red for inhibitory
    1 1 1;  % White for no connection
    0 0 1   % Blue for excitatory
];

figure
imagesc(connMatrixPlot, [1 3]);
colormap(cmap);
colorbar('Ticks', [1, 2, 3], ...
         'TickLabels', {'Inhibitory', 'No Connection', 'Excitatory'});
xlabel('Presynaptic Neuron Index');
ylabel('Postsynaptic Neuron Index');
title('Ground-Truth Cross-TT Connectivity Matrix');
axis image

%catch 
%disp('Error processing i number. Writing cells as not classified')
%end

%% Save Review Session
pause(1)

firing_rates = All_Cells_combined(iii).classific_firingRate;

review_pairs = initialize_review_pairs(exc_pairs_all, inh_pairs_all, rebound_override_pairs, spikes_full.shankID);

review_session = struct();
review_session.version = 2;
review_session.session_index = iii;
review_session.session_label = sprintf('Animal %s Day %s Session %s', ...
    num2str(sessInfo(iii).animal), num2str(sessInfo(iii).day), num2str(iii));
review_session.mainDir = sessInfo(iii).mainDir;
review_session.animal = sessInfo(iii).animal;
review_session.day = sessInfo(iii).day;
review_session.behavior_session_dirs = sessInfo(iii).sessDirs(1:number_OFs);
review_session.n_units = numel(mono_res_excitatory.n);
review_session.tt_files = tt_files;
review_session.firing_rates = firing_rates;
review_session.auto_cell_classification = cell_classification;
review_session.behavior_spike_times = behavior_spike_times;
review_session.exc_pairs_all = exc_pairs_all;
review_session.inh_pairs_all = inh_pairs_all;
review_session.rebound_override_pairs = rebound_override_pairs;
review_session.review_pairs = review_pairs;
review_session.spike_times = allOF_tSp;
review_session.mono_res_excitatory = mono_res_excitatory;
review_session.mono_res_inhibitory = mono_res_inhibitory;
review_session.plot_config = struct( ...
    'sig_window_exc_ms', sig_window_exc_ms, ...
    'sig_window_inh_ms', sig_window_inh_ms, ...
    'plot_bin_full_s', plot_bin_full_s, ...
    'plot_bin_zoom_exc_s', plot_bin_zoom_exc_s, ...
    'plot_bin_zoom_inh_s', plot_bin_zoom_inh_s, ...
    'plot_bin_raw_s', plot_bin_raw_s, ...
    'zoom_window_exc_ms', zoom_window_exc_ms, ...
    'zoom_window_inh_ms', zoom_window_inh_ms, ...
    'full_window_ms', [min(ccg_ts_ms) max(ccg_ts_ms)] ...
    );
review_session.detection_config = detection_config;
review_session.settings_signature = settings_signature;
review_session.auto_pair_counts = struct( ...
    'total_connections_tested', total_pairs_evaluated, ...
    'total_pairs_evaluated', total_pairs_evaluated, ...
    'total_pairs_crossTT', total_pairs_crossTT, ...
    'total_pairs_sameTT', total_pairs_sameTT, ...
    'identified_candidate_pairs', numel(review_pairs), ...
    'detected_exc_all', size(exc_pairs_all,1), ...
    'detected_inh_all', size(inh_pairs_all,1), ...
    'detected_exc_crossTT', size(exc_pairs_crossTT,1), ...
    'detected_exc_sameTT', size(exc_pairs_sameTT,1), ...
    'detected_inh_crossTT', size(inh_pairs_crossTT,1), ...
    'detected_inh_sameTT', size(inh_pairs_sameTT,1) ...
    );

review_session = merge_review_session_state(review_session, review_file);
review_session.classifier_run_mode = review_generation_mode;
save(review_file, 'review_session', '-v7.3');
fprintf('Saved review session to %s\n', review_file);
processed_session_count = processed_session_count + 1;

cell_classification_plot = cell_classification;

%% GUI for manual curation
%gui_MonoSyn(mono_res_excitatory) % Shows the GUI for manual curation

%% Plot some strongly connected cross-correlograms
if make_plots ==1
outfolder = ('W:\Haseeb\Pictures\CellPairs');
if ~exist(outfolder, 'dir')
    mkdir(outfolder)
end
outfolder_sameTT = fullfile(outfolder, 'SameTT');
outfolder_crossTT = fullfile(outfolder, 'CrossTT');
outfolder_rebound = fullfile(outfolder, 'ReboundOverrides');
if ~exist(outfolder_sameTT, 'dir')
    mkdir(outfolder_sameTT)
end
if ~exist(outfolder_crossTT, 'dir')
    mkdir(outfolder_crossTT)
end
if ~exist(outfolder_rebound, 'dir')
    mkdir(outfolder_rebound)
end


connections = [];

if ~isempty(accepted_exc_pairs_crossTT)
    connections = [connections; accepted_exc_pairs_crossTT, repmat(1, size(accepted_exc_pairs_crossTT,1), 1), zeros(size(accepted_exc_pairs_crossTT,1), 1)];
end

if ~isempty(accepted_exc_pairs_sameTT)
    connections = [connections; accepted_exc_pairs_sameTT, repmat(1, size(accepted_exc_pairs_sameTT,1), 1), ones(size(accepted_exc_pairs_sameTT,1), 1)];
end

if ~isempty(accepted_inh_pairs_crossTT)
    connections = [connections; accepted_inh_pairs_crossTT, repmat(2, size(accepted_inh_pairs_crossTT,1), 1), zeros(size(accepted_inh_pairs_crossTT,1), 1)];
end

if ~isempty(accepted_inh_pairs_sameTT)
    connections = [connections; accepted_inh_pairs_sameTT, repmat(2, size(accepted_inh_pairs_sameTT,1), 1), ones(size(accepted_inh_pairs_sameTT,1), 1)];
end


for idx = 1:size(connections,1)
    presyn = connections(idx,1);
    postsyn = connections(idx,2);
    connection_type = connections(idx,3);
    same_tt_pair = connections(idx,4) == 1;

    if connection_type == 3
        continue
    end

    if ~(cell_classification_plot(presyn) == 1 || cell_classification_plot(presyn) == 2)
        continue
    end

     if connection_type == 2
        mono_plot = mono_res_inhibitory;
        sig_window_this_ms = sig_window_inh_ms;
        zoom_window_this_ms = zoom_window_inh_ms;
        plot_bin_zoom_this_s = plot_bin_zoom_inh_s;
    else
        mono_plot = mono_res_excitatory;
        sig_window_this_ms = sig_window_exc_ms;
        zoom_window_this_ms = zoom_window_exc_ms;
        plot_bin_zoom_this_s = plot_bin_zoom_exc_s;
    end
    zoom_lag_this_s = zoom_window_this_ms / 1000;

    plot_method_suffix = 'CellExplorerMatched';
    % Fine-grid predictor/bounds from CellExplorer result (keep detection output as is)
    pred_plot_fine = squeeze(mono_plot.Pred(:,presyn,postsyn));
    upper_plot_fine = squeeze(mono_plot.Bounds(:,presyn,postsyn,1));
    lower_plot_fine = squeeze(mono_plot.Bounds(:,presyn,postsyn,2));
    pred_ts_ms_fine = mono_plot.ccgTs(:) * 1000;

    % -------- FULL PLOT (rebinned detector CCG) --------
    full_range_ms = [min(pred_ts_ms_fine) max(pred_ts_ms_fine)];

    [ccg_full_ts_ms, ccg_full_values] = rebin_ccg_counts(ccg_ts_ms(:), ...
        mono_plot.ccgR(:,presyn,postsyn), plot_bin_full_s * 1000, full_range_ms);
    plot_method_label_full = sprintf('Detector CCG (rebinned to %.1f ms)', plot_bin_full_s * 1000);

    [pred_full_ts_ms, pred_full] = rebin_ccg_counts(pred_ts_ms_fine, pred_plot_fine, plot_bin_full_s * 1000, full_range_ms);
    [~, upper_full] = rebin_ccg_counts(pred_ts_ms_fine, upper_plot_fine, plot_bin_full_s * 1000, full_range_ms);
    [~, lower_full] = rebin_ccg_counts(pred_ts_ms_fine, lower_plot_fine, plot_bin_full_s * 1000, full_range_ms);

    connection_label = {'excitatory', 'inhibitory', 'both'};
    presyn_tt = tt_files{presyn};
    postsyn_tt = tt_files{postsyn};

    if same_tt_pair
        pair_scope_text = 'same-TT';
        pair_scope_suffix = 'sameTT';
        outfolder_this = outfolder_sameTT;
    else
        pair_scope_text = 'cross-TT';
        pair_scope_suffix = 'crossTT';
        outfolder_this = outfolder_crossTT;
    end

    plot_ccg_panel(ccg_full_ts_ms, ccg_full_values, ...
        pred_full_ts_ms, pred_full, upper_full, lower_full, ...
        full_range_ms, sig_window_this_ms, presyn_tt, postsyn_tt, firing_rates, presyn, postsyn, ...
        connection_label{connection_type}, plot_method_label_full, sprintf('(full, %s pair)', pair_scope_text));

    filename = fullfile(outfolder_this, sprintf('Animal%s-day%s-Cell-%s-%s-%s-%s-%s.pdf', ...
        num2str(sessInfo(iii).animal), num2str(sessInfo(iii).day), ...
        presyn_tt(1:end-2), postsyn_tt(1:end-2), ...
        connection_label{connection_type}, pair_scope_suffix, plot_method_suffix));
    exportgraphics(gcf, filename, 'ContentType', 'vector');
    filename_png = strrep(filename, '.pdf', '.png');
    exportgraphics(gcf, filename_png, 'Resolution', 300);

    % -------- ZOOM PLOT (connection-specific display bins) --------
    [ccg_zoom_ts_ms, ccg_zoom_values] = rebin_ccg_counts(ccg_ts_ms(:), ...
        mono_plot.ccgR(:,presyn,postsyn), plot_bin_zoom_this_s*1000, zoom_window_this_ms);
    plot_method_label_zoom = sprintf('Detector CCG (%.1f ms bins)', plot_bin_zoom_this_s * 1000);

    [pred_zoom_ts_ms, pred_zoom] = rebin_ccg_counts(pred_ts_ms_fine, pred_plot_fine, ...
        plot_bin_zoom_this_s*1000, zoom_window_this_ms);
    [~, upper_zoom] = rebin_ccg_counts(pred_ts_ms_fine, upper_plot_fine, ...
        plot_bin_zoom_this_s*1000, zoom_window_this_ms);
    [~, lower_zoom] = rebin_ccg_counts(pred_ts_ms_fine, lower_plot_fine, ...
        plot_bin_zoom_this_s*1000, zoom_window_this_ms);

    plot_ccg_panel(ccg_zoom_ts_ms, ccg_zoom_values, ...
        pred_zoom_ts_ms, pred_zoom, upper_zoom, lower_zoom, ...
        zoom_window_this_ms, sig_window_this_ms, presyn_tt, postsyn_tt, firing_rates, presyn, postsyn, ...
        connection_label{connection_type}, plot_method_label_zoom, ...
        sprintf('(zoom %g to %+g ms, %s pair)', zoom_window_this_ms(1), zoom_window_this_ms(2), pair_scope_text));

    filename_zoom = fullfile(outfolder_this, sprintf('Animal%s-day%s-Cell-%s-%s-%s-%s-%s-zoom.pdf', ...
        num2str(sessInfo(iii).animal), num2str(sessInfo(iii).day), ...
        presyn_tt(1:end-2), postsyn_tt(1:end-2), ...
        connection_label{connection_type}, pair_scope_suffix, plot_method_suffix));
    exportgraphics(gcf, filename_zoom, 'ContentType', 'vector');
    filename_zoom_png = strrep(filename_zoom, '.pdf', '.png');
    exportgraphics(gcf, filename_zoom_png, 'Resolution', 300);

    % -------- EXTRA RAW 0.1 MS PLOT (display-only) --------
    [ccg_raw_values, ccg_raw_ts] = CrossCorrel(allOF_tSp{presyn}, allOF_tSp{postsyn}, ...
        plot_bin_raw_s, zoom_lag_this_s);
    ccg_raw_values = double(ccg_raw_values(:));
    ccg_raw_ts_ms = ccg_raw_ts(:) * 1000;

    plot_ccg_panel(ccg_raw_ts_ms, ccg_raw_values, ...
        [], [], [], [], ...
        zoom_window_this_ms, sig_window_this_ms, presyn_tt, postsyn_tt, firing_rates, presyn, postsyn, ...
        connection_label{connection_type}, 'Raw CCG (0.1 ms, display-only)', ...
        sprintf('(zoom %g to %+g ms, %s pair)', zoom_window_this_ms(1), zoom_window_this_ms(2), pair_scope_text));

    filename_raw = fullfile(outfolder_this, sprintf('Animal%s-day%s-Cell-%s-%s-%s-%s-Raw0100us.pdf', ...
        num2str(sessInfo(iii).animal), num2str(sessInfo(iii).day), ...
        presyn_tt(1:end-2), postsyn_tt(1:end-2), ...
        connection_label{connection_type}, pair_scope_suffix));
    exportgraphics(gcf, filename_raw, 'ContentType', 'vector');
    filename_raw_png = strrep(filename_raw, '.pdf', '.png');
    exportgraphics(gcf, filename_raw_png, 'Resolution', 300);

    % -------- EXTRA RESIDUAL PLOT FOR INHIBITORY PAIRS --------
    if connection_type == 2 && save_inhibitory_residual_plots
        residual_zoom = ccg_zoom_values - pred_zoom;
        upper_resid_zoom = upper_zoom - pred_zoom;
        lower_resid_zoom = lower_zoom - pred_zoom;

        plot_ccg_residual_panel(ccg_zoom_ts_ms, residual_zoom, upper_resid_zoom, lower_resid_zoom, ...
            zoom_window_this_ms, sig_window_this_ms, presyn_tt, postsyn_tt, firing_rates, presyn, postsyn, ...
            connection_label{connection_type}, plot_method_label_zoom, sprintf('(residual: CCG - Pred, %s pair)', pair_scope_text));

        filename_resid = fullfile(outfolder_this, sprintf('Animal%s-day%s-Cell-%s-%s-%s-%s-%s-residual.pdf', ...
            num2str(sessInfo(iii).animal), num2str(sessInfo(iii).day), ...
            presyn_tt(1:end-2), postsyn_tt(1:end-2), ...
            connection_label{connection_type}, pair_scope_suffix, plot_method_suffix));
        exportgraphics(gcf, filename_resid, 'ContentType', 'vector');
        filename_resid_png = strrep(filename_resid, '.pdf', '.png');
        exportgraphics(gcf, filename_resid_png, 'Resolution', 300);
    end

    close all
    
    
end

for idx = 1:size(rebound_override_pairs,1)
    presyn = rebound_override_pairs(idx,1);
    postsyn = rebound_override_pairs(idx,2);

    same_tt_pair = spikes_full.shankID(presyn) == spikes_full.shankID(postsyn);
    if same_tt_pair
        pair_scope_text = 'same-TT';
        pair_scope_suffix = 'sameTT';
    else
        pair_scope_text = 'cross-TT';
        pair_scope_suffix = 'crossTT';
    end

    mono_plot = mono_res_inhibitory;
    zoom_window_this_ms = zoom_window_inh_ms;
    sig_window_this_ms = sig_window_inh_ms;
    plot_bin_zoom_this_s = plot_bin_zoom_inh_s;

    pred_plot_fine = squeeze(mono_plot.Pred(:,presyn,postsyn));
    upper_plot_fine = squeeze(mono_plot.Bounds(:,presyn,postsyn,1));
    lower_plot_fine = squeeze(mono_plot.Bounds(:,presyn,postsyn,2));
    pred_ts_ms_fine = mono_plot.ccgTs(:) * 1000;

    [first_trough_ms, later_peak_ms] = get_rebound_timing(mono_plot, presyn, postsyn);

    [ccg_zoom_ts_ms, ccg_zoom_values] = rebin_ccg_counts(ccg_ts_ms(:), ...
        mono_plot.ccgR(:,presyn,postsyn), plot_bin_zoom_this_s*1000, zoom_window_this_ms);
    [pred_zoom_ts_ms, pred_zoom] = rebin_ccg_counts(pred_ts_ms_fine, pred_plot_fine, ...
        plot_bin_zoom_this_s*1000, zoom_window_this_ms);
    [~, upper_zoom] = rebin_ccg_counts(pred_ts_ms_fine, upper_plot_fine, ...
        plot_bin_zoom_this_s*1000, zoom_window_this_ms);
    [~, lower_zoom] = rebin_ccg_counts(pred_ts_ms_fine, lower_plot_fine, ...
        plot_bin_zoom_this_s*1000, zoom_window_this_ms);

    presyn_tt = tt_files{presyn};
    postsyn_tt = tt_files{postsyn};

    title_suffix = sprintf('(reassigned excit->inh, %s pair, first trough %.1f ms, later peak %.1f ms)', ...
        pair_scope_text, first_trough_ms, later_peak_ms);

    plot_ccg_panel(ccg_zoom_ts_ms, ccg_zoom_values, ...
        pred_zoom_ts_ms, pred_zoom, upper_zoom, lower_zoom, ...
        zoom_window_this_ms, sig_window_this_ms, presyn_tt, postsyn_tt, firing_rates, presyn, postsyn, ...
        'inhibitory rebound override', sprintf('Detector CCG (%.1f ms bins)', plot_bin_zoom_this_s * 1000), title_suffix);

    filename_rebound = fullfile(outfolder_rebound, sprintf('Animal%s-day%s-Cell-%s-%s-%s-ReboundOverride-Detector.pdf', ...
        num2str(sessInfo(iii).animal), num2str(sessInfo(iii).day), ...
        presyn_tt(1:end-2), postsyn_tt(1:end-2), pair_scope_suffix));
    exportgraphics(gcf, filename_rebound, 'ContentType', 'vector');
    filename_rebound_png = strrep(filename_rebound, '.pdf', '.png');
    exportgraphics(gcf, filename_rebound_png, 'Resolution', 300);

    [ccg_raw_values, ccg_raw_ts] = CrossCorrel(allOF_tSp{presyn}, allOF_tSp{postsyn}, ...
        plot_bin_raw_s, zoom_window_this_ms / 1000);
    ccg_raw_values = double(ccg_raw_values(:));
    ccg_raw_ts_ms = ccg_raw_ts(:) * 1000;

    plot_ccg_panel(ccg_raw_ts_ms, ccg_raw_values, ...
        [], [], [], [], ...
        zoom_window_this_ms, sig_window_this_ms, presyn_tt, postsyn_tt, firing_rates, presyn, postsyn, ...
        'inhibitory rebound override', 'Raw CCG (0.1 ms, display-only)', title_suffix);

    filename_rebound_raw = fullfile(outfolder_rebound, sprintf('Animal%s-day%s-Cell-%s-%s-%s-ReboundOverride-Raw0100us.pdf', ...
        num2str(sessInfo(iii).animal), num2str(sessInfo(iii).day), ...
        presyn_tt(1:end-2), postsyn_tt(1:end-2), pair_scope_suffix));
    exportgraphics(gcf, filename_rebound_raw, 'ContentType', 'vector');
    filename_rebound_raw_png = strrep(filename_rebound_raw, '.pdf', '.png');
    exportgraphics(gcf, filename_rebound_raw_png, 'Resolution', 300);

    residual_zoom = ccg_zoom_values - pred_zoom;
    upper_resid_zoom = upper_zoom - pred_zoom;
    lower_resid_zoom = lower_zoom - pred_zoom;

    plot_ccg_residual_panel(ccg_zoom_ts_ms, residual_zoom, upper_resid_zoom, lower_resid_zoom, ...
        zoom_window_this_ms, sig_window_this_ms, presyn_tt, postsyn_tt, firing_rates, presyn, postsyn, ...
        'inhibitory rebound override', sprintf('Detector CCG (%.1f ms bins)', plot_bin_zoom_this_s * 1000), title_suffix);

    filename_rebound_resid = fullfile(outfolder_rebound, sprintf('Animal%s-day%s-Cell-%s-%s-%s-ReboundOverride-Residual.pdf', ...
        num2str(sessInfo(iii).animal), num2str(sessInfo(iii).day), ...
        presyn_tt(1:end-2), postsyn_tt(1:end-2), pair_scope_suffix));
    exportgraphics(gcf, filename_rebound_resid, 'ContentType', 'vector');
    filename_rebound_resid_png = strrep(filename_rebound_resid, '.pdf', '.png');
    exportgraphics(gcf, filename_rebound_resid_png, 'Resolution', 300);

    close all
end

clear cell_classification exc_cells inh_cells_all cell_classification_plot
end


end

fprintf('CCG classifier finished. Processed sessions: %d | Skipped existing review files: %d | Mode: %s\n', ...
    processed_session_count, skipped_existing_review_count, review_generation_mode);
end

function config = load_classification_config(config_path)

if nargin < 1 || isempty(config_path)
    config_path = fullfile(fileparts(mfilename('fullpath')), 'classification_config.json');
end

config = struct('cellExplorerPath', '', 'mclustPath', '');
if ~exist(config_path, 'file')
    return
end

try
    decoded = jsondecode(fileread(config_path));
catch ME
    warning('Could not read classification config %s: %s', config_path, ME.message);
    return
end

if isfield(decoded, 'cellExplorerPath')
    config.cellExplorerPath = char(string(decoded.cellExplorerPath));
end
if isfield(decoded, 'mclustPath')
    config.mclustPath = char(string(decoded.mclustPath));
end
end

function add_repo_matlab_code_path()

matlabCodeFolder = fullfile(fileparts(mfilename('fullpath')), '..');
if exist(matlabCodeFolder, 'dir')
    addpath(genpath(matlabCodeFolder));
end
end

function add_configured_dependency_paths(config, custom_settings)

cellExplorerPath = char(string(get_override_value(custom_settings, 'cellExplorerPath', config.cellExplorerPath)));
mclustPath = char(string(get_override_value(custom_settings, 'mclustPath', config.mclustPath)));

add_dependency_path(cellExplorerPath, 'CellExplorer');
add_dependency_path(mclustPath, 'MClust');

additionalPaths = get_override_value(custom_settings, 'additionalPaths', {});
if ischar(additionalPaths) || isstring(additionalPaths)
    additionalPaths = cellstr(string(additionalPaths));
end
for iPath = 1:numel(additionalPaths)
    add_dependency_path(additionalPaths{iPath}, 'additional dependency');
end
end

function add_dependency_path(path_value, dependency_name)

path_value = char(string(path_value));
if isempty(path_value)
    return
end

if exist(path_value, 'dir')
    addpath(genpath(path_value));
else
    warning('%s path does not exist: %s', dependency_name, path_value);
end
end

function sessionInfoPath = resolve_session_info_path()

repoRoot = fullfile(fileparts(mfilename('fullpath')), '..', '..');
candidatePaths = { ...
    fullfile(repoRoot, 'sessionInfo.mat'), ...
    fullfile(repoRoot, 'Data', 'sessionInfo.mat'), ...
    fullfile(repoRoot, 'Analysis_scripts', 'DataOrganization', 'sessionInfo.mat'), ...
    fullfile(pwd, 'sessionInfo.mat')};

sessionInfoPath = first_existing_file(candidatePaths);
if isempty(sessionInfoPath)
    error('Could not find sessionInfo.mat. Pass custom_settings.sessionInfoPath or place it in the repository root/Data folder.')
end
end

function allCellsPath = resolve_all_cells_path()

repoRoot = fullfile(fileparts(mfilename('fullpath')), '..', '..');
candidatePaths = { ...
    fullfile(repoRoot, 'All_Cells_combined.mat'), ...
    fullfile(repoRoot, 'Data', 'All_Cells_combined.mat'), ...
    fullfile(pwd, 'All_Cells_combined.mat')};

allCellsPath = first_existing_file(candidatePaths);
if isempty(allCellsPath)
    error('Could not find All_Cells_combined.mat. Pass custom_settings.allCellsPath or place it in the repository root/Data folder.')
end
end

function reviewRoot = default_review_root()

repoRoot = fullfile(fileparts(mfilename('fullpath')), '..', '..');
reviewRoot = fullfile(repoRoot, 'Data', 'CCG_review_sessions');
end

function pathOut = first_existing_file(candidatePaths)

pathOut = '';
for iPath = 1:numel(candidatePaths)
    if exist(candidatePaths{iPath}, 'file')
        pathOut = candidatePaths{iPath};
        return
    end
end
end

function value = get_override_value(settings_struct, field_name, default_value)

value = default_value;
if nargin < 1 || isempty(settings_struct) || ~isstruct(settings_struct)
    return
end

if isfield(settings_struct, field_name) && ~isempty(settings_struct.(field_name))
    value = settings_struct.(field_name);
end
end

function session_indices = default_ccg_session_indices()

session_indices = [11:17, 21:25, 31:38, 41:47, 52:56, 81:82, 84:88, 101:110, 121, ...
    141:153, 161:171, 181:199, 201:219, 222:239, 241:256];
end

function review_generation_mode = resolve_review_generation_mode(default_mode)

review_generation_mode = default_mode;
if nargin < 1 || isempty(review_generation_mode)
    review_generation_mode = 'only_missing_review_files';
end

if ~usejava('desktop')
    return
end

selection_choice = questdlg( ...
    ['Run all selected sessions again, or only create review files for sessions ' ...
     'that do not have one yet?'], ...
    'CCG Review File Mode', ...
    'Only Missing Review Files', 'Re-test All Selected Sessions', 'Cancel', ...
    'Only Missing Review Files');

if isempty(selection_choice) || strcmp(selection_choice, 'Cancel')
    review_generation_mode = '';
elseif strcmp(selection_choice, 'Re-test All Selected Sessions')
    review_generation_mode = 'retest_all_selected_sessions';
else
    review_generation_mode = 'only_missing_review_files';
end
end

function tf = should_skip_review_session(review_file, review_generation_mode)

tf = strcmp(review_generation_mode, 'only_missing_review_files') && exist(review_file, 'file') == 2;
end

function [cor, lags] = CrossCorrel(ts1, ts2, binsize, lag)
ac = isequal(ts1,ts2);
ts1 = ts1(:);
ts2 = ts2(:);
if ~issorted(ts1)
    ts1 = sort(ts1);
end
if ~issorted(ts2)
    ts2 = sort(ts2);
end

lags = lag(1)+binsize/2:binsize:lag(2)+binsize/2;
lags = lags(:);
cor = zeros(size(lags));
if isempty(ts1) || isempty(ts2)
    return
end

db = nan(length(ts1), 3);
psth = [];
s1 = 1;
spkind = 1;
while spkind <= length(ts1)
    s = s1;
    while s < length(ts2) && ts2(s) < ts1(spkind)+lag(1)
        s = s+1;
    end
    s1 = s;
    f = s;
    while f < length(ts2) && ts2(f) <= ts1(spkind)+lag(2)
        f = f+1;
    end
    if ts2(s)<=ts1(spkind)+lag(2)
        db(spkind, :) = [s f-1 ts1(spkind)];
    end
    spkind = spkind+1;
end

valid_rows = ~isnan(db(:,1));
if ~any(valid_rows)
    return
end
db_valid = db(valid_rows,:);
dspk = diff(db_valid(:,1:2), 1, 2);
if isempty(dspk)
    return
end
for i = 0:max(dspk)
    where = dspk>=i;
    tf = db_valid(:,1);
    tf = tf(where)+i;
    psth = [psth; ts2(tf)-db_valid(where,3)];
end
if ac
    psth(psth==0) = [];
end
cor = hist(psth, lags);
cor = cor(:);
end

function [t_out_ms, y_out] = rebin_ccg_counts(t_in_ms, y_in, new_bin_ms, t_range_ms)
% Rebin counts/predictor/bounds using only integer multiples of original bins.

t_in_ms = t_in_ms(:);
y_in = y_in(:);

old_bin_ms = median(diff(t_in_ms));
ratio = new_bin_ms / old_bin_ms;

assert(abs(ratio - round(ratio)) < 1e-9, ...
    'Display bin size must be an integer multiple of original bin size.');

ratio = round(ratio);

keep = t_in_ms >= t_range_ms(1) & t_in_ms <= t_range_ms(2);
t = t_in_ms(keep);
y = y_in(keep);

n_keep = floor(numel(y)/ratio) * ratio;
t = t(1:n_keep);
y = y(1:n_keep);

y_out = sum(reshape(y, ratio, []), 1)';
t_out_ms = mean(reshape(t, ratio, []), 1)';
end


function plot_ccg_panel(bar_t_ms, bar_y, pred_t_ms, pred_y, upper_y, lower_y, ...
    xwin_ms, sig_window_ms, presyn_tt, postsyn_tt, firing_rates, presyn, postsyn, ...
    connection_label_text, plot_method_label, title_suffix)

figure
bar(bar_t_ms, bar_y, 1, 'FaceColor', [0.3 0.3 0.3], 'EdgeColor', 'none')
hold on

yl = ylim;
h_sig = patch([sig_window_ms(1) sig_window_ms(2) sig_window_ms(2) sig_window_ms(1)], ...
              [yl(1) yl(1) yl(2) yl(2)], ...
              [0.2 0.8 0.2], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
uistack(h_sig, 'bottom');

if ~isempty(pred_t_ms)
    plot(pred_t_ms, pred_y, 'k-', 'LineWidth', 1.2)
    plot(pred_t_ms, upper_y, 'k--', 'LineWidth', 0.8)
    plot(pred_t_ms, lower_y, 'k--', 'LineWidth', 0.8)
end

xlim(xwin_ms)
xlabel('time [ms]');
ylabel('Count');

title(sprintf(['Pre: %s (FR: %s Hz) - Post: %s (FR: %s Hz)\n' ...
    'Connection type: %s | %s %s'], ...
    replace(presyn_tt(1:end-2), '_', ' '), num2str(firing_rates(presyn)), ...
    replace(postsyn_tt(1:end-2), '_', ' '), num2str(firing_rates(postsyn)), ...
    connection_label_text, plot_method_label, title_suffix));

hold off
end

function plot_ccg_residual_panel(bar_t_ms, residual_y, upper_resid, lower_resid, ...
    xwin_ms, sig_window_ms, presyn_tt, postsyn_tt, firing_rates, presyn, postsyn, ...
    connection_label_text, plot_method_label, title_suffix)

figure
bar(bar_t_ms, residual_y, 1, 'FaceColor', [0.3 0.3 0.3], 'EdgeColor', 'none')
hold on

yl = ylim;
h_sig = patch([sig_window_ms(1) sig_window_ms(2) sig_window_ms(2) sig_window_ms(1)], ...
              [yl(1) yl(1) yl(2) yl(2)], ...
              [0.2 0.8 0.2], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
uistack(h_sig, 'bottom');

plot(bar_t_ms, zeros(size(bar_t_ms)), 'k-', 'LineWidth', 1.0)
plot(bar_t_ms, upper_resid, 'k--', 'LineWidth', 0.8)
plot(bar_t_ms, lower_resid, 'k--', 'LineWidth', 0.8)

xlim(xwin_ms)
xlabel('time [ms]');
ylabel('CCG - Pred');

title(sprintf(['Pre: %s (FR: %s Hz) - Post: %s (FR: %s Hz)\n' ...
    'Connection type: %s | %s %s'], ...
    replace(presyn_tt(1:end-2), '_', ' '), num2str(firing_rates(presyn)), ...
    replace(postsyn_tt(1:end-2), '_', ' '), num2str(firing_rates(postsyn)), ...
    connection_label_text, plot_method_label, title_suffix));

hold off
end

function mono_res = expand_mono_res_to_full_indexing(mono_res, valid_units, nUnits_full)

orig_ids = find(valid_units);
nValid = numel(orig_ids);

% Store mapping explicitly
mono_res.valid_units = valid_units;
mono_res.orig_ids = orig_ids;

% Expand n back to full length
if isfield(mono_res, 'n') && ~isempty(mono_res.n)
    n_full = zeros(nUnits_full,1);
    n_full(orig_ids) = mono_res.n(:);
    mono_res.n = n_full;
else
    mono_res.n = zeros(nUnits_full,1);
end

% Remap significant connection indices back to original unit IDs
if isfield(mono_res, 'sig_con_excitatory')
    mono_res.sig_con_excitatory = normalize_connection_pairs(mono_res.sig_con_excitatory);
    if ~isempty(mono_res.sig_con_excitatory)
        mono_res.sig_con_excitatory = orig_ids(mono_res.sig_con_excitatory);
        mono_res.sig_con_excitatory = normalize_connection_pairs(mono_res.sig_con_excitatory);
    end
else
    mono_res.sig_con_excitatory = zeros(0,2);
end

if isfield(mono_res, 'sig_con_inhibitory')
    mono_res.sig_con_inhibitory = normalize_connection_pairs(mono_res.sig_con_inhibitory);
    if ~isempty(mono_res.sig_con_inhibitory)
        mono_res.sig_con_inhibitory = orig_ids(mono_res.sig_con_inhibitory);
        mono_res.sig_con_inhibitory = normalize_connection_pairs(mono_res.sig_con_inhibitory);
    end
else
    mono_res.sig_con_inhibitory = zeros(0,2);
end

% Expand ccgR back to full [lags x units x units]
if isfield(mono_res, 'ccgR') && ~isempty(mono_res.ccgR)
    sz = size(mono_res.ccgR);
    ccgR_full = zeros(sz(1), nUnits_full, nUnits_full, 'like', mono_res.ccgR);
    ccgR_full(:, orig_ids, orig_ids) = mono_res.ccgR;
    mono_res.ccgR = ccgR_full;
end

% Expand Pred back to full [lags x units x units]
if isfield(mono_res, 'Pred') && ~isempty(mono_res.Pred)
    sz = size(mono_res.Pred);
    Pred_full = zeros(sz(1), nUnits_full, nUnits_full, 'like', mono_res.Pred);
    Pred_full(:, orig_ids, orig_ids) = mono_res.Pred;
    mono_res.Pred = Pred_full;
end

% Expand Bounds back to full [lags x units x units x 2]
if isfield(mono_res, 'Bounds') && ~isempty(mono_res.Bounds)
    sz = size(mono_res.Bounds);
    Bounds_full = zeros(sz(1), nUnits_full, nUnits_full, sz(4), 'like', mono_res.Bounds);
    Bounds_full(:, orig_ids, orig_ids, :) = mono_res.Bounds;
    mono_res.Bounds = Bounds_full;
end
end

function conn = normalize_connection_pairs(conn)
% Ensure connection list is always N x 2
% Accepts:
%   []        -> []
%   N x 2     -> unchanged
%   2 x N     -> transposed to N x 2
%   2 x 1     -> reshaped to 1 x 2
%   1 x 2     -> unchanged

if isempty(conn)
    conn = zeros(0,2);
    return
end

sz = size(conn);

if isequal(sz, [1 2]) || (numel(sz) == 2 && sz(2) == 2)
    % already N x 2
    return
elseif isequal(sz, [2 1])
    conn = reshape(conn, 1, 2);
elseif sz(1) == 2 && sz(2) ~= 2
    conn = conn.';
elseif numel(conn) == 2
    conn = reshape(conn, 1, 2);
else
    error('Connection array has unexpected size: [%s]', num2str(size(conn)));
end
end

function [crossTT_pairs, sameTT_pairs] = split_pairs_by_tetrode(pairs, shank_ids)

if isempty(pairs)
    crossTT_pairs = zeros(0,2);
    sameTT_pairs = zeros(0,2);
    return
end

same_tt_mask = shank_ids(pairs(:,1)) == shank_ids(pairs(:,2));
sameTT_pairs = pairs(same_tt_mask, :);
crossTT_pairs = pairs(~same_tt_mask, :);
end

function [total_pairs_evaluated, total_pairs_crossTT, total_pairs_sameTT] = ...
        compute_total_pair_counts(shank_ids, valid_units)

valid_idx = find(valid_units(:));
if numel(valid_idx) < 2
    total_pairs_evaluated = 0;
    total_pairs_crossTT = 0;
    total_pairs_sameTT = 0;
    return
end

valid_shanks = shank_ids(valid_idx);
same_mask = valid_shanks(:) == valid_shanks(:)';
same_mask(1:size(same_mask,1)+1:end) = false;

total_pairs_sameTT = sum(same_mask(:));
total_pairs_evaluated = numel(valid_idx) * (numel(valid_idx) - 1);
total_pairs_crossTT = total_pairs_evaluated - total_pairs_sameTT;
end

function rebound_pairs = detect_inhibitory_rebound_pairs(mono_res, exc_pairs)

if isempty(exc_pairs)
    rebound_pairs = zeros(0,2);
    return
end

causal_mask = mono_res.ccgTs(:) >= 0.001 & mono_res.ccgTs(:) <= 0.004;
rebound_pairs = zeros(0,2);

for idx = 1:size(exc_pairs,1)
    pre = exc_pairs(idx,1);
    post = exc_pairs(idx,2);

    cch = squeeze(mono_res.ccgR(:,pre,post));
    upper = squeeze(mono_res.Bounds(:,pre,post,1));
    lower = squeeze(mono_res.Bounds(:,pre,post,2));

    cch_causal = cch(causal_mask);
    upper_causal = upper(causal_mask);
    lower_causal = lower(causal_mask);

    sig_pos = cch_causal > upper_causal;
    sig_neg = cch_causal < lower_causal;
    sig_any = sig_pos | sig_neg;

    if ~any(sig_any)
        continue
    end

    first_sig_idx = find(sig_any, 1, 'first');

    if sig_neg(first_sig_idx) && any(sig_pos(first_sig_idx+1:end))
        rebound_pairs = [rebound_pairs; pre post];
    end
end

rebound_pairs = unique_pair_rows(rebound_pairs);
end

function [first_trough_ms, later_peak_ms] = get_rebound_timing(mono_res, pre, post)

first_trough_ms = nan;
later_peak_ms = nan;

causal_mask = mono_res.ccgTs(:) >= 0.001 & mono_res.ccgTs(:) <= 0.004;
causal_ts_ms = mono_res.ccgTs(causal_mask) * 1000;

cch = squeeze(mono_res.ccgR(:,pre,post));
upper = squeeze(mono_res.Bounds(:,pre,post,1));
lower = squeeze(mono_res.Bounds(:,pre,post,2));

cch_causal = cch(causal_mask);
upper_causal = upper(causal_mask);
lower_causal = lower(causal_mask);

sig_pos = cch_causal > upper_causal;
sig_neg = cch_causal < lower_causal;
sig_any = sig_pos | sig_neg;

if ~any(sig_any)
    return
end

first_sig_idx = find(sig_any, 1, 'first');
if ~sig_neg(first_sig_idx)
    return
end

first_trough_ms = causal_ts_ms(first_sig_idx);
later_peak_rel_idx = find(sig_pos(first_sig_idx+1:end), 1, 'first');

if isempty(later_peak_rel_idx)
    return
end

later_peak_ms = causal_ts_ms(first_sig_idx + later_peak_rel_idx);
end

function pairs_out = remove_pair_rows(pairs_in, pairs_to_remove)

if isempty(pairs_in)
    pairs_out = zeros(0,2);
    return
end

if isempty(pairs_to_remove)
    pairs_out = pairs_in;
    return
end

remove_mask = ismember(pairs_in, pairs_to_remove, 'rows');
pairs_out = pairs_in(~remove_mask, :);
end

function pairs_out = unique_pair_rows(pairs_in)

if isempty(pairs_in)
    pairs_out = zeros(0,2);
    return
end

pairs_out = unique(pairs_in, 'rows', 'stable');
end

function review_session = merge_review_session_state(review_session, review_file)

review_session.review_pairs = standardize_review_pairs(review_session.review_pairs);
review_session.version = 2;
review_session.candidate_signature = make_candidate_signature(review_session.review_pairs);
review_session.total_connections_tested = numel(review_session.review_pairs);
review_session.total_connections_verified = sum([review_session.review_pairs.verified]);
review_session.review_complete = review_session.total_connections_verified == review_session.total_connections_tested;
review_session.writeback_complete = false;
review_session.last_writeback = '';
review_session.curated_file = '';
review_session.reviewed_cell_classification = [];
review_session.previous_writeback_complete = false;
review_session.previous_last_writeback = '';
review_session.previous_curated_file = '';
review_session.last_classifier_run = char(datetime('now'));
review_session.merge_summary = struct( ...
    'matched_pairs', 0, ...
    'preserved_verified_pairs', 0, ...
    'reset_due_to_candidate_change', 0, ...
    'new_pairs', review_session.total_connections_tested, ...
    'settings_matched_previous', false ...
    );

if ~exist(review_file, 'file')
    return
end

loaded_existing = load(review_file, 'review_session');
if ~isfield(loaded_existing, 'review_session')
    return
end

existing_review = standardize_review_session(loaded_existing.review_session);
same_settings = isfield(existing_review, 'settings_signature') && ...
    strcmp(existing_review.settings_signature, review_session.settings_signature);

if isfield(existing_review, 'writeback_complete')
    review_session.previous_writeback_complete = existing_review.writeback_complete;
end
if isfield(existing_review, 'last_writeback')
    review_session.previous_last_writeback = existing_review.last_writeback;
end
if isfield(existing_review, 'curated_file')
    review_session.previous_curated_file = existing_review.curated_file;
end

review_session.merge_summary.settings_matched_previous = same_settings;

if same_settings
    [review_session.review_pairs, merge_summary] = merge_review_pairs(review_session.review_pairs, existing_review.review_pairs);
    review_session.merge_summary = merge_summary;
    review_session.merge_summary.settings_matched_previous = true;
else
    review_session.review_pairs = reset_review_pairs(review_session.review_pairs);
    review_session.merge_summary.reset_due_to_candidate_change = review_session.total_connections_tested;
end

review_session.review_pairs = standardize_review_pairs(review_session.review_pairs);
review_session.candidate_signature = make_candidate_signature(review_session.review_pairs);
review_session.total_connections_tested = numel(review_session.review_pairs);
review_session.total_connections_verified = sum([review_session.review_pairs.verified]);
review_session.review_complete = review_session.total_connections_verified == review_session.total_connections_tested;

keep_writeback = same_settings && ...
    isfield(existing_review, 'candidate_signature') && ...
    strcmp(existing_review.candidate_signature, review_session.candidate_signature) && ...
    isfield(existing_review, 'writeback_complete') && existing_review.writeback_complete && ...
    review_session.review_complete;

if keep_writeback
    review_session.writeback_complete = true;
    if isfield(existing_review, 'last_writeback')
        review_session.last_writeback = existing_review.last_writeback;
    end
    if isfield(existing_review, 'curated_file')
        review_session.curated_file = existing_review.curated_file;
    end
    if isfield(existing_review, 'reviewed_cell_classification')
        review_session.reviewed_cell_classification = existing_review.reviewed_cell_classification;
    end
end
end

function review_session = standardize_review_session(review_session)

if ~isfield(review_session, 'review_pairs') || isempty(review_session.review_pairs)
    if ~isfield(review_session, 'n_units') || isempty(review_session.n_units)
        shank_stub = 1;
    else
        shank_stub = zeros(max(review_session.n_units, 1), 1);
    end
    review_session.review_pairs = initialize_review_pairs(zeros(0,2), zeros(0,2), zeros(0,2), shank_stub);
end

review_session.review_pairs = standardize_review_pairs(review_session.review_pairs);

if ~isfield(review_session, 'candidate_signature')
    review_session.candidate_signature = make_candidate_signature(review_session.review_pairs);
end
if ~isfield(review_session, 'settings_signature')
    review_session.settings_signature = '';
end
if ~isfield(review_session, 'total_connections_tested')
    review_session.total_connections_tested = numel(review_session.review_pairs);
end
if ~isfield(review_session, 'total_connections_verified')
    review_session.total_connections_verified = sum([review_session.review_pairs.verified]);
end
if ~isfield(review_session, 'review_complete')
    review_session.review_complete = review_session.total_connections_verified == review_session.total_connections_tested;
end
if ~isfield(review_session, 'writeback_complete')
    review_session.writeback_complete = isfield(review_session, 'last_writeback') && ~isempty(review_session.last_writeback);
end
if ~isfield(review_session, 'last_writeback')
    review_session.last_writeback = '';
end
if ~isfield(review_session, 'curated_file')
    review_session.curated_file = '';
end
if ~isfield(review_session, 'reviewed_cell_classification')
    review_session.reviewed_cell_classification = [];
end
if ~isfield(review_session, 'previous_writeback_complete')
    review_session.previous_writeback_complete = false;
end
if ~isfield(review_session, 'previous_last_writeback')
    review_session.previous_last_writeback = '';
end
if ~isfield(review_session, 'previous_curated_file')
    review_session.previous_curated_file = '';
end
end

function review_pairs = standardize_review_pairs(review_pairs)

if isempty(review_pairs)
    review_pairs = initialize_review_pairs(zeros(0,2), zeros(0,2), zeros(0,2), 1);
    return
end

for idx = 1:numel(review_pairs)
    if ~isfield(review_pairs, 'verified') || isempty(review_pairs(idx).verified)
        review_pairs(idx).verified = false;
    else
        review_pairs(idx).verified = logical(review_pairs(idx).verified);
    end

    if ~isfield(review_pairs, 'last_reviewed') || isempty(review_pairs(idx).last_reviewed)
        review_pairs(idx).last_reviewed = '';
    end

    if ~isfield(review_pairs, 'current_decision') || isempty(review_pairs(idx).current_decision)
        review_pairs(idx).current_decision = review_pairs(idx).initial_decision;
    end
end
end

function [merged_pairs, merge_summary] = merge_review_pairs(new_pairs, old_pairs)

new_pairs = standardize_review_pairs(new_pairs);
old_pairs = standardize_review_pairs(old_pairs);

merge_summary = struct( ...
    'matched_pairs', 0, ...
    'preserved_verified_pairs', 0, ...
    'reset_due_to_candidate_change', 0, ...
    'new_pairs', 0, ...
    'settings_matched_previous', true ...
    );

if isempty(new_pairs)
    merged_pairs = new_pairs;
    return
end

merged_pairs = new_pairs;

for idx = 1:numel(new_pairs)
    match_idx = find(arrayfun(@(p) p.pre == new_pairs(idx).pre && p.post == new_pairs(idx).post, old_pairs), 1, 'first');

    if isempty(match_idx)
        merge_summary.new_pairs = merge_summary.new_pairs + 1;
        continue
    end

    merge_summary.matched_pairs = merge_summary.matched_pairs + 1;

    old_pair = old_pairs(match_idx);
    merged_pairs(idx).current_decision = old_pair.current_decision;
    merged_pairs(idx).last_reviewed = old_pair.last_reviewed;

    same_candidate = old_pair.auto_type == new_pairs(idx).auto_type && ...
        old_pair.same_tt_pair == new_pairs(idx).same_tt_pair && ...
        old_pair.rebound_override == new_pairs(idx).rebound_override;

    if same_candidate
        merged_pairs(idx).verified = old_pair.verified;
        if old_pair.verified
            merge_summary.preserved_verified_pairs = merge_summary.preserved_verified_pairs + 1;
        end
    else
        merged_pairs(idx).verified = false;
        merge_summary.reset_due_to_candidate_change = merge_summary.reset_due_to_candidate_change + 1;
    end
end
end

function review_pairs = reset_review_pairs(review_pairs)

review_pairs = standardize_review_pairs(review_pairs);

for idx = 1:numel(review_pairs)
    review_pairs(idx).current_decision = review_pairs(idx).initial_decision;
    review_pairs(idx).verified = false;
    review_pairs(idx).last_reviewed = '';
end
end

function signature = make_settings_signature(detection_config)

signature = sprintf(['bin=%.7f|duration=%.7f|excSig=%.7f|inhSig=%.7f|alpha=%.7f|' ...
    'conv=%.7f|ref=%.7f,%.7f|causalStart=%.7f|mask=%.7f,%.7f|rebound=%.7f,%.7f|sr=%d'], ...
    detection_config.bin_s, detection_config.duration_s, ...
    detection_config.sig_window_exc_s, detection_config.sig_window_inh_s, ...
    detection_config.alpha, detection_config.conv_w_s, ...
    detection_config.reference_window_s(1), detection_config.reference_window_s(2), ...
    detection_config.causal_window_start_s, ...
    detection_config.same_shank_mask_s(1), detection_config.same_shank_mask_s(2), ...
    detection_config.rebound_window_s(1), detection_config.rebound_window_s(2), ...
    detection_config.sr);
end

function signature = make_candidate_signature(review_pairs)

review_pairs = standardize_review_pairs(review_pairs);

if isempty(review_pairs)
    signature = 'none';
    return
end

tokens = cell(numel(review_pairs), 1);
for idx = 1:numel(review_pairs)
    tokens{idx} = sprintf('%d-%d-%d-%d-%d', ...
        review_pairs(idx).pre, review_pairs(idx).post, review_pairs(idx).auto_type, ...
        review_pairs(idx).same_tt_pair, review_pairs(idx).rebound_override);
end

signature = strjoin(tokens, ';');
end

function review_pairs = initialize_review_pairs(exc_pairs, inh_pairs, rebound_pairs, shank_ids)

all_pairs = unique_pair_rows([exc_pairs; inh_pairs]);

empty_entry = struct( ...
    'pre', 0, ...
    'post', 0, ...
    'auto_type', 0, ...
    'initial_decision', 0, ...
    'current_decision', 0, ...
    'same_tt_pair', false, ...
    'rebound_override', false, ...
    'verified', false, ...
    'last_reviewed', '' ...
    );

if isempty(all_pairs)
    review_pairs = repmat(empty_entry, 0, 1);
    return
end

review_pairs = repmat(empty_entry, size(all_pairs,1), 1);

for idx = 1:size(all_pairs,1)
    pre = all_pairs(idx,1);
    post = all_pairs(idx,2);

    is_exc = ismember([pre post], exc_pairs, 'rows');
    is_inh = ismember([pre post], inh_pairs, 'rows');
    is_rebound = ismember([pre post], rebound_pairs, 'rows');

    if is_exc && is_inh
        auto_type = 3;
        initial_decision = 0;
    elseif is_exc
        auto_type = 1;
        initial_decision = 1;
    elseif is_inh
        auto_type = 2;
        initial_decision = 2;
    else
        auto_type = 0;
        initial_decision = 0;
    end

    review_pairs(idx).pre = pre;
    review_pairs(idx).post = post;
    review_pairs(idx).auto_type = auto_type;
    review_pairs(idx).initial_decision = initial_decision;
    review_pairs(idx).current_decision = initial_decision;
    review_pairs(idx).same_tt_pair = shank_ids(pre) == shank_ids(post);
    review_pairs(idx).rebound_override = is_rebound;
    review_pairs(idx).verified = false;
    review_pairs(idx).last_reviewed = '';
end
end
