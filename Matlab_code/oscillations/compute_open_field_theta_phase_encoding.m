function compute_open_field_theta_phase_encoding(varargin)
%COMPUTE_OPEN_FIELD_THETA_PHASE_ENCODING Compute theta phase locking during OF encoding.
%
% Purpose
%   Assign the instantaneous theta-band LFP phase to open-field spikes and
%   summarize each cell's phase locking across the selected OF sessions.
%   The function does not define cFos/engram membership; downstream group
%   comparisons should read those labels from All_Cells_combined.optotagged.
%
% Inputs
%   sessionInfoPath  .mat file containing sessInfo. If empty, the function
%                    resolves the standard sessionInfo.mat location.
%   allCellsPath     .mat file containing All_Cells_combined. If empty, the
%                    script resolves the standard All_Cells_combined.mat.
%   sessions         Numeric session indices to process. Empty processes all
%                    sessions present in both sessInfo and All_Cells_combined.
%   openFieldFolders OF folders to pool, default {'of1','of2','of3'}.
%   thetaBand        LFP bandpass range in Hz, default [6 14].
%   speedThreshold   Minimum tracking speed for included spikes, default 2 cm/s.
%   minSpikes        Minimum included spikes per cell for circular metrics.
%
% Required session data
%   sessInfo(iii).mainDir, .tList, and .cellLayerChann.
%   processedData/indata_of.mat with tracking time and speed.
%   OF spike files readable by readSpikeDataOnly/fixSpikes.
%   OF EEG/LFP files readable by readCRTsd/Data/Range.
%
% Outputs
%   All_Cells_combined(iii).thetaPhase_pref     preferred phase in radians.
%   All_Cells_combined(iii).thetaPhase_R        mean resultant length.
%   All_Cells_combined(iii).thetaPhase_p        Rayleigh p-value.
%   All_Cells_combined(iii).thetaPhase_z        Rayleigh z statistic.
%   All_Cells_combined(iii).thetaPhase_nSpikes  included spike count.
%
%   Per-session details are saved to processedData/OF_theta_phase.mat,
%   including OF-specific phase vectors, spike times, metrics, and metadata.
%   All_Cells_combined is updated on disk when saveAllCells is true.

p = inputParser;
p.addParameter('configPath', '', @(x) ischar(x) || isstring(x));
p.addParameter('sessionInfoPath', '', @(x) ischar(x) || isstring(x));
p.addParameter('allCellsPath', '', @(x) ischar(x) || isstring(x));
p.addParameter('mclustPath', '', @(x) ischar(x) || isstring(x));
p.addParameter('additionalPaths', {});
p.addParameter('sessions', [], @(x) isempty(x) || isnumeric(x));
p.addParameter('openFieldFolders', {'of1','of2','of3'});
p.addParameter('thetaBand', [6 14], @(x) isnumeric(x) && numel(x) == 2);
p.addParameter('speedThreshold', 2, @(x) isnumeric(x) && isscalar(x));
p.addParameter('minSpikes', 5, @(x) isnumeric(x) && isscalar(x));
p.addParameter('downsamplingFactor', 16, @(x) isnumeric(x) && isscalar(x));
p.addParameter('minThetaAmplitudePercentile', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
p.addParameter('saveAllCells', true, @(x) islogical(x) || isnumeric(x));
p.parse(varargin{:});
opts = p.Results;

opts.configPath = char(string(opts.configPath));
opts.mclustPath = char(string(opts.mclustPath));
opts.additionalPaths = normalize_path_list(opts.additionalPaths);
opts.openFieldFolders = cellstr(string(opts.openFieldFolders));
opts.thetaBand = double(opts.thetaBand(:)');
opts.speedThreshold = double(opts.speedThreshold);
opts.minSpikes = double(opts.minSpikes);
opts.downsamplingFactor = double(opts.downsamplingFactor);
opts.saveAllCells = logical(opts.saveAllCells);

config = load_oscillation_config(opts.configPath);
add_dependency_paths(config, opts);

sessionInfoPath = char(opts.sessionInfoPath);
if isempty(sessionInfoPath)
    sessionInfoPath = resolve_session_info_path();
end

allCellsPath = char(opts.allCellsPath);
if isempty(allCellsPath)
    allCellsPath = resolve_all_cells_path();
end

loadedSessionInfo = load(sessionInfoPath, 'sessInfo');
sessInfo = loadedSessionInfo.sessInfo;
loadedAllCells = load(allCellsPath, 'All_Cells_combined');
All_Cells_combined = loadedAllCells.All_Cells_combined;

if isempty(opts.sessions)
    sessionsToRun = 1:min(numel(sessInfo), numel(All_Cells_combined));
else
    sessionsToRun = opts.sessions(:)';
end

fprintf('\n=== OF theta phase encoding ===\n');
fprintf('SessionInfo: %s\n', sessionInfoPath);
fprintf('All_Cells_combined: %s\n', allCellsPath);
fprintf('Theta band: %.1f-%.1f Hz | speed > %.2f cm/s | min spikes: %d\n', ...
    opts.thetaBand(1), opts.thetaBand(2), opts.speedThreshold, opts.minSpikes);

for iSession = 1:numel(sessionsToRun)
    iii = sessionsToRun(iSession);
    if iii < 1 || iii > numel(sessInfo) || iii > numel(All_Cells_combined)
        fprintf('Skipping session index %d: outside available session range.\n', iii);
        continue
    end

    if ~isfield(sessInfo(iii), 'mainDir') || ~isfield(sessInfo(iii), 'tList') || ...
            ~isfield(sessInfo(iii), 'cellLayerChann')
        fprintf('Skipping session %d: missing mainDir, tList, or cellLayerChann.\n', iii);
        continue
    end

    mainDir = normalize_text_field(sessInfo(iii).mainDir);
    tListName = normalize_text_field(sessInfo(iii).tList);
    if isempty(mainDir) || isempty(tListName)
        fprintf('Skipping session %d: empty mainDir or tList.\n', iii);
        continue
    end

    tListPath = fullfile(mainDir, tListName);
    if ~exist(tListPath, 'file')
        fprintf('Skipping session %d: missing tList file %s.\n', iii, tListPath);
        continue
    end

    indataFile = fullfile(mainDir, 'processedData', 'indata_of.mat');
    if ~exist(indataFile, 'file')
        fprintf('Skipping session %d: missing OF tracking file %s.\n', iii, indataFile);
        continue
    end

    ttFiles = read_tt_list(tListPath);
    if isempty(ttFiles)
        fprintf('Skipping session %d: tList contains no cells.\n', iii);
        continue
    end

    nCells = max(numel(ttFiles), max_num_cells_in_session(All_Cells_combined(iii)));
    if nCells == 0
        fprintf('Skipping session %d: no cells found.\n', iii);
        continue
    end

    loadedIndata = load(indataFile, 'indata');
    indataAll = loadedIndata.indata;
    pooledPhasesByCell = cell(nCells, 1);
    ofSessions = repmat(empty_of_result(nCells), numel(opts.openFieldFolders), 1);

    fprintf('\nSession %d/%d (i=%d, animal=%s, day=%s): %d cells\n', ...
        iSession, numel(sessionsToRun), iii, ...
        normalize_text_field(get_struct_field(sessInfo(iii), 'animal')), ...
        normalize_text_field(get_struct_field(sessInfo(iii), 'day')), nCells);

    for iOF = 1:numel(opts.openFieldFolders)
        ofFolder = opts.openFieldFolders{iOF};
        ofNumber = parse_of_number(ofFolder, iOF);

        if ofNumber < 1 || ofNumber > numel(indataAll)
            ofSessions(iOF) = empty_of_result(nCells, ofFolder, 'skipped', 'Missing matching indata entry.');
            fprintf('  %s skipped: missing matching indata entry.\n', ofFolder);
            continue
        end

        if ~exist(fullfile(mainDir, ofFolder), 'dir')
            ofSessions(iOF) = empty_of_result(nCells, ofFolder, 'skipped', 'Missing OF folder.');
            fprintf('  %s skipped: missing OF folder.\n', ofFolder);
            continue
        end

        try
            ofSessions(iOF) = compute_open_field_phase_session( ...
                sessInfo(iii), mainDir, ofFolder, indataAll(ofNumber), ttFiles, nCells, opts);
        catch ME
            ofSessions(iOF) = empty_of_result(nCells, ofFolder, 'error', ME.message);
            fprintf('  %s error: %s\n', ofFolder, ME.message);
            continue
        end

        for c = 1:nCells
            pooledPhasesByCell{c} = [pooledPhasesByCell{c}; ofSessions(iOF).phaseByCell{c}(:)];
        end

        fprintf('  %s done: %d running spikes assigned to theta phase.\n', ...
            ofFolder, sum(ofSessions(iOF).nSpikes(isfinite(ofSessions(iOF).nSpikes))));
    end

    [thetaPhase_pref, thetaPhase_R, thetaPhase_p, thetaPhase_z, thetaPhase_nSpikes] = ...
        compute_phase_metrics_by_cell(pooledPhasesByCell, opts.minSpikes);

    All_Cells_combined(iii).thetaPhase_pref = thetaPhase_pref;
    All_Cells_combined(iii).thetaPhase_R = thetaPhase_R;
    All_Cells_combined(iii).thetaPhase_p = thetaPhase_p;
    All_Cells_combined(iii).thetaPhase_z = thetaPhase_z;
    All_Cells_combined(iii).thetaPhase_nSpikes = thetaPhase_nSpikes;

    phase_metadata = build_phase_metadata(sessInfo(iii), iii, ttFiles, opts);
    of_theta_phase = struct();
    of_theta_phase.metadata = phase_metadata;
    of_theta_phase.ofSessions = ofSessions;
    of_theta_phase.summary.thetaPhase_pref = thetaPhase_pref;
    of_theta_phase.summary.thetaPhase_R = thetaPhase_R;
    of_theta_phase.summary.thetaPhase_p = thetaPhase_p;
    of_theta_phase.summary.thetaPhase_z = thetaPhase_z;
    of_theta_phase.summary.thetaPhase_nSpikes = thetaPhase_nSpikes;

    processedDir = fullfile(mainDir, 'processedData');
    if ~exist(processedDir, 'dir')
        mkdir(processedDir);
    end
    save(fullfile(processedDir, 'OF_theta_phase.mat'), 'of_theta_phase', ...
        'thetaPhase_pref', 'thetaPhase_R', 'thetaPhase_p', 'thetaPhase_z', 'thetaPhase_nSpikes');

    fprintf('  Saved compact summary to All_Cells_combined and details to %s.\n', ...
        fullfile(processedDir, 'OF_theta_phase.mat'));
end

if opts.saveAllCells
    save(allCellsPath, 'All_Cells_combined', '-append');
    fprintf('\nSaved updated All_Cells_combined: %s\n', allCellsPath);
end

fprintf('OF theta phase encoding completed.\n');
end


function ofResult = compute_open_field_phase_session(sessEntry, mainDir, ofFolder, indata, ttFiles, nCells, opts)

lfp = load_eeg_local(mainDir, ofFolder, sessEntry.cellLayerChann, opts.downsamplingFactor);
[lfp.ts, lfp.samp] = align_lfp_to_tracking(indata.t, lfp.ts, lfp.samp);

thetaSignal = filter_theta_lfp(lfp.samp, lfp.sampFreq, opts.thetaBand);
thetaSignal = zscore_vector(thetaSignal);
analyticSignal = hilbert(thetaSignal);
thetaPhase = angle(analyticSignal);
thetaAmplitude = abs(analyticSignal);

speedAtLfp = interp1(indata.t(:), indata.v(:), lfp.ts(:), 'linear', NaN);
runningLfpMask = speedAtLfp > opts.speedThreshold;

thetaAmplitudeThreshold = NaN;
if ~isempty(opts.minThetaAmplitudePercentile)
    thetaAmplitudeThreshold = percentile_local(thetaAmplitude(runningLfpMask), opts.minThetaAmplitudePercentile);
end

spikeData = readSpikeDataOnly(fullfile(mainDir, ofFolder), ttFiles);
spikeTimes = fixSpikes(spikeData);

phaseByCell = cell(nCells, 1);
spikeTimesByCell = cell(nCells, 1);
for c = 1:nCells
    phaseByCell{c} = [];
    spikeTimesByCell{c} = [];
end

for c = 1:min(numel(spikeTimes), nCells)
    cellSpikes = spikeTimes{c};
    if isempty(cellSpikes)
        continue
    end

    cellSpikes = double(cellSpikes(:));
    spikeSpeed = interp1(indata.t(:), indata.v(:), cellSpikes, 'linear', NaN);
    validSpikes = cellSpikes >= lfp.ts(1) & cellSpikes <= lfp.ts(end) & ...
        spikeSpeed > opts.speedThreshold;

    if ~isempty(opts.minThetaAmplitudePercentile)
        spikeThetaAmplitude = interp1(lfp.ts(:), thetaAmplitude(:), cellSpikes, 'linear', NaN);
        validSpikes = validSpikes & spikeThetaAmplitude >= thetaAmplitudeThreshold;
    end

    validSpikeTimes = cellSpikes(validSpikes);
    if isempty(validSpikeTimes)
        continue
    end

    spikePhases = interp1(lfp.ts(:), thetaPhase(:), validSpikeTimes, 'nearest', NaN);
    validPhase = isfinite(spikePhases);
    phaseByCell{c} = spikePhases(validPhase);
    phaseByCell{c} = phaseByCell{c}(:);
    spikeTimesByCell{c} = validSpikeTimes(validPhase);
    spikeTimesByCell{c} = spikeTimesByCell{c}(:);
end

[pref, R, p, z, nSpikes] = compute_phase_metrics_by_cell(phaseByCell, opts.minSpikes);

ofResult = empty_of_result(nCells, ofFolder, 'ok', '');
ofResult.thetaBand = opts.thetaBand;
ofResult.speedThreshold = opts.speedThreshold;
ofResult.minSpikes = opts.minSpikes;
ofResult.minThetaAmplitudePercentile = opts.minThetaAmplitudePercentile;
ofResult.thetaAmplitudeThreshold = thetaAmplitudeThreshold;
ofResult.lfpSampFreq = lfp.sampFreq;
ofResult.nValidLfpSamples = nnz(runningLfpMask);
ofResult.phaseByCell = phaseByCell;
ofResult.spikeTimesByCell = spikeTimesByCell;
ofResult.pref = pref;
ofResult.R = R;
ofResult.p = p;
ofResult.z = z;
ofResult.nSpikes = nSpikes;
end


function ofResult = empty_of_result(nCells, folder, status, message)

if nargin < 1 || isempty(nCells)
    nCells = 0;
end
if nargin < 2
    folder = '';
end
if nargin < 3
    status = '';
end
if nargin < 4
    message = '';
end

ofResult = struct();
ofResult.folder = char(folder);
ofResult.status = char(status);
ofResult.message = char(message);
ofResult.thetaBand = [NaN NaN];
ofResult.speedThreshold = NaN;
ofResult.minSpikes = NaN;
ofResult.minThetaAmplitudePercentile = [];
ofResult.thetaAmplitudeThreshold = NaN;
ofResult.lfpSampFreq = NaN;
ofResult.nValidLfpSamples = NaN;
ofResult.phaseByCell = cell(nCells, 1);
ofResult.spikeTimesByCell = cell(nCells, 1);
ofResult.pref = nan(nCells, 1);
ofResult.R = nan(nCells, 1);
ofResult.p = nan(nCells, 1);
ofResult.z = nan(nCells, 1);
ofResult.nSpikes = nan(nCells, 1);
end


function [pref, R, p, z, nSpikes] = compute_phase_metrics_by_cell(phaseByCell, minSpikes)

nCells = numel(phaseByCell);
pref = nan(nCells, 1);
R = nan(nCells, 1);
p = nan(nCells, 1);
z = nan(nCells, 1);
nSpikes = zeros(nCells, 1);

for c = 1:nCells
    phases = phaseByCell{c};
    phases = phases(:);
    phases = phases(isfinite(phases));
    n = numel(phases);
    nSpikes(c) = n;

    if n < minSpikes
        continue
    end

    sumSin = sum(sin(phases));
    sumCos = sum(cos(phases));
    pref(c) = atan2(sumSin, sumCos);
    R(c) = sqrt(sumSin.^2 + sumCos.^2) / n;
    z(c) = n * R(c).^2;
    p(c) = exp(sqrt(1 + 4*n + 4*(n.^2 - (n*R(c)).^2)) - (1 + 2*n));
end
end


function lfp = load_eeg_local(infolder, task, channel, downsamplingFactor)

if iscell(channel)
    channel = channel{1};
end
if ischar(channel) || isstring(channel)
    channel = str2double(char(string(channel)));
end
channel = double(channel(1));
eegfilename = fullfile(infolder, task, strcat('CSC', num2str(channel), '.ncs'));
fprintf('Reading: %s\n', eegfilename);
[eeg, sFreq] = readCRTsd(eegfilename);

lfp.samp = Data(eeg) * -1;
lfp.samp = downsample(lfp.samp(:), downsamplingFactor);
lfp.ts = Range(eeg);
lfp.ts = downsample(lfp.ts(:), downsamplingFactor);
lfp.sampFreq = sFreq / downsamplingFactor;
end


function [eeg_ts, eeg_raw] = align_lfp_to_tracking(trackingTime, eeg_ts, eeg_raw)

if eeg_ts(end) / trackingTime(end) > 10
    eeg_ts = eeg_ts * 1e-4;
end

[~, startIdx] = min(abs(trackingTime(1) - eeg_ts));
[~, stopIdx] = min(abs(trackingTime(end) - eeg_ts));

eeg_raw = eeg_raw(startIdx:stopIdx);
eeg_ts = eeg_ts(startIdx:stopIdx);
end


function thetaSignal = filter_theta_lfp(lfpSignal, sampFreq, thetaBand)

lfpSignal = double(lfpSignal(:));
badSamples = ~isfinite(lfpSignal);
if any(badSamples)
    goodIdx = find(~badSamples);
    if isempty(goodIdx)
        error('LFP signal contains no finite samples.');
    end
    lfpSignal(badSamples) = interp1(goodIdx, lfpSignal(goodIdx), find(badSamples), 'linear', 'extrap');
end

try
    thetaSignal = thetaphase.BandpassFilter(lfpSignal, sampFreq, thetaBand);
catch
    [b, a] = butter(3, thetaBand ./ (sampFreq / 2), 'bandpass');
    thetaSignal = filtfilt(b, a, lfpSignal);
end

thetaSignal = thetaSignal(:);
end


function x = zscore_vector(x)

x = double(x(:));
mu = mean(x(isfinite(x)));
sigma = std(x(isfinite(x)));
if ~isfinite(sigma) || sigma == 0
    sigma = 1;
end
x = (x - mu) ./ sigma;
end


function value = percentile_local(x, prctileValue)

x = sort(x(isfinite(x)));
if isempty(x)
    value = NaN;
    return
end

prctileValue = max(0, min(100, prctileValue));
idx = 1 + (numel(x) - 1) * prctileValue / 100;
lo = floor(idx);
hi = ceil(idx);
if lo == hi
    value = x(lo);
else
    value = x(lo) + (idx - lo) * (x(hi) - x(lo));
end
end


function phase_metadata = build_phase_metadata(sessEntry, sessionIndex, ttFiles, opts)

phase_metadata = struct();
phase_metadata.session_index = sessionIndex;
phase_metadata.animal = get_struct_field(sessEntry, 'animal');
phase_metadata.day = get_struct_field(sessEntry, 'day');
phase_metadata.mainDir = normalize_text_field(get_struct_field(sessEntry, 'mainDir'));
phase_metadata.t_list = ttFiles(:);
phase_metadata.t_list_file = normalize_text_field(get_struct_field(sessEntry, 'tList'));
phase_metadata.thetaBand = opts.thetaBand;
phase_metadata.speedThreshold = opts.speedThreshold;
phase_metadata.minSpikes = opts.minSpikes;
phase_metadata.downsamplingFactor = opts.downsamplingFactor;
phase_metadata.minThetaAmplitudePercentile = opts.minThetaAmplitudePercentile;
phase_metadata.openFieldFolders = opts.openFieldFolders;
phase_metadata.phaseConvention = 'angle(hilbert(zscored 6-14 Hz LFP)); LFP sign follows Data*-1 loading convention';
phase_metadata.source_script = mfilename;
phase_metadata.generated_on = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
end


function ttFiles = read_tt_list(tListPath)

fid = fopen(tListPath);
if fid < 0
    error('Could not open tList file: %s', tListPath);
end

ttFiles = {};
while true
    tline = fgetl(fid);
    if ~ischar(tline)
        break
    end
    tline = strtrim(tline);
    if ~isempty(tline)
        ttFiles{end+1, 1} = tline; %#ok<AGROW>
    end
end
fclose(fid);
end


function nCells = max_num_cells_in_session(S)

nCells = 0;
candidateFields = {'thetaPhase_pref', 'thetaMod_score', 'optotagged', 'SD_0_ms_delay', ...
    'GMM_based_classification_days', 'final_classification_numeric', 'of_avg_fir_rate'};

for iField = 1:numel(candidateFields)
    fieldName = candidateFields{iField};
    if isfield(S, fieldName) && ~isempty(S.(fieldName))
        nCells = max(nCells, numel(S.(fieldName)));
    end
end
end


function ofNumber = parse_of_number(ofFolder, fallback)

tokens = regexp(char(ofFolder), 'of(\d+)', 'tokens', 'once');
if isempty(tokens)
    ofNumber = fallback;
else
    ofNumber = str2double(tokens{1});
end
end


function txt = normalize_text_field(value)

if nargin == 0 || isempty(value)
    txt = '';
    return
end

if iscell(value)
    value = value{1};
end

if isstring(value)
    txt = char(value(1));
elseif ischar(value)
    txt = value;
elseif isnumeric(value) && isscalar(value)
    txt = num2str(value);
else
    txt = char(string(value));
end
end


function value = get_struct_field(S, fieldName)

if isfield(S, fieldName)
    value = S.(fieldName);
else
    value = [];
end
end


function sessionInfoPath = resolve_session_info_path()

repoRoot = fullfile(fileparts(mfilename('fullpath')), '..', '..');
candidates = { ...
    fullfile(repoRoot, 'sessionInfo.mat'), ...
    fullfile(repoRoot, 'Data', 'sessionInfo.mat'), ...
    fullfile(repoRoot, 'Analysis_scripts', 'DataOrganization', 'sessionInfo.mat'), ...
    fullfile(pwd, 'sessionInfo.mat'), ...
    fullfile(pwd, 'DataOrganization', 'sessionInfo.mat')};
sessionInfoPath = resolve_existing_file(candidates, 'sessionInfo.mat');
end


function allCellsPath = resolve_all_cells_path()

repoRoot = fullfile(fileparts(mfilename('fullpath')), '..', '..');
candidates = { ...
    fullfile(repoRoot, 'All_Cells_combined.mat'), ...
    fullfile(repoRoot, 'Data', 'All_Cells_combined.mat'), ...
    fullfile(pwd, 'All_Cells_combined.mat')};
allCellsPath = resolve_existing_file(candidates, 'All_Cells_combined.mat');
end


function filePath = resolve_existing_file(candidates, description)

filePath = '';
for iCandidate = 1:numel(candidates)
    candidate = candidates{iCandidate};
    if exist(candidate, 'file')
        filePath = candidate;
        return
    end
end

error('Could not find %s. Pass an explicit path as a name-value argument.', description);
end


function config = load_oscillation_config(configPath)

if nargin < 1 || isempty(configPath)
    configPath = fullfile(fileparts(mfilename('fullpath')), '..', 'classification', 'classification_config.json');
end

config = struct('mclustPath', '', 'thetaAnalysisPath', '', 'oscillationAnalysisPath', '');
if exist(configPath, 'file') ~= 2
    return
end

try
    decoded = jsondecode(fileread(configPath));
catch ME
    warning('Could not read config %s: %s', configPath, ME.message);
    return
end

if isfield(decoded, 'mclustPath')
    config.mclustPath = char(string(decoded.mclustPath));
end
if isfield(decoded, 'thetaAnalysisPath')
    config.thetaAnalysisPath = char(string(decoded.thetaAnalysisPath));
end
if isfield(decoded, 'oscillationAnalysisPath')
    config.oscillationAnalysisPath = char(string(decoded.oscillationAnalysisPath));
end
end


function add_dependency_paths(config, opts)

matlabCodeFolder = fullfile(fileparts(mfilename('fullpath')), '..');
if exist(matlabCodeFolder, 'dir') == 7
    addpath(genpath(matlabCodeFolder));
end

mclustPath = char(string(get_option_value(opts, 'mclustPath', config.mclustPath)));
add_dependency_path(mclustPath, 'MClust');
add_dependency_path(config.thetaAnalysisPath, 'theta analysis dependency');
add_dependency_path(config.oscillationAnalysisPath, 'oscillation analysis dependency');

for iPath = 1:numel(opts.additionalPaths)
    add_dependency_path(opts.additionalPaths{iPath}, 'additional dependency');
end
end


function add_dependency_path(folderPath, dependencyName)

folderPath = char(string(folderPath));
if isempty(folderPath)
    return
end

if exist(folderPath, 'dir') == 7
    addpath(genpath(folderPath));
else
    warning('%s path does not exist: %s', dependencyName, folderPath);
end
end


function values = normalize_path_list(values)

if ischar(values) || isstring(values)
    values = cellstr(string(values));
end
end


function value = get_option_value(settingsStruct, fieldName, defaultValue)

value = defaultValue;
if nargin < 1 || isempty(settingsStruct) || ~isstruct(settingsStruct)
    return
end

if isfield(settingsStruct, fieldName) && ~isempty(settingsStruct.(fieldName))
    value = settingsStruct.(fieldName);
end
end
