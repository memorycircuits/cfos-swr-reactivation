function rebuild_All_Cells_combined_CCG_classification_from_curated(all_cells_file, curated_session_root)

if nargin < 1 || isempty(all_cells_file)
    all_cells_file = resolve_all_cells_path();
end

if nargin < 2 || isempty(curated_session_root)
    curated_session_root = resolve_curated_session_root();
end

loaded = load(all_cells_file, 'All_Cells_combined');
if ~isfield(loaded, 'All_Cells_combined')
    error('All_Cells_combined was not found in %s.', all_cells_file);
end

All_Cells_combined = loaded.All_Cells_combined;
n_sessions = numel(All_Cells_combined);

for idx = 1:n_sessions
    unit_count = infer_all_cells_session_unit_count(All_Cells_combined(idx));
    All_Cells_combined(idx).CCGbased_CellClassfication = zeros(unit_count, 1);
end

listing = dir(fullfile(curated_session_root, '*-Curated.mat'));
if isempty(listing)
    warning('No curated session files were found in %s. All sessions were set to zero labels only.', curated_session_root);
    save(all_cells_file, 'All_Cells_combined', '-v7.3');
    return
end

updated_sessions = [];
length_adjusted_sessions = [];
skipped_sessions = [];

for idx = 1:numel(listing)
    loaded_curated = load(fullfile(listing(idx).folder, listing(idx).name), 'curated_session');
    if ~isfield(loaded_curated, 'curated_session')
        skipped_sessions(end+1,1) = idx; %#ok<AGROW>
        continue
    end

    curated_session = loaded_curated.curated_session;
    session_index = curated_session.session_index;

    if session_index < 1 || session_index > n_sessions
        skipped_sessions(end+1,1) = session_index; %#ok<AGROW>
        continue
    end

    base_count = infer_all_cells_session_unit_count(All_Cells_combined(session_index));
    curated_count = infer_curated_session_unit_count(curated_session);
    final_count = max([base_count, curated_count, numel(curated_session.final_cell_classification)]);

    if final_count == 0
        All_Cells_combined(session_index).CCGbased_CellClassfication = zeros(0, 1);
        updated_sessions(end+1,1) = session_index; %#ok<AGROW>
        continue
    end

    classification = zeros(final_count, 1);
    source_classification = curated_session.final_cell_classification(:);

    if ~isempty(source_classification)
        classification(1:numel(source_classification)) = source_classification;
    end

    if base_count ~= final_count || curated_count ~= numel(source_classification)
        length_adjusted_sessions(end+1,1) = session_index; %#ok<AGROW>
    end

    All_Cells_combined(session_index).CCGbased_CellClassfication = classification;
    updated_sessions(end+1,1) = session_index; %#ok<AGROW>
end

save(all_cells_file, 'All_Cells_combined', '-v7.3');

fprintf('Rebuilt All_Cells_combined.CCGbased_CellClassfication from curated sessions.\n');
fprintf('Source curated folder: %s\n', curated_session_root);
fprintf('Sessions updated from curated files: %d\n', numel(updated_sessions));
fprintf('Sessions with zero-only labels (no curated file): %d\n', n_sessions - numel(unique(updated_sessions)));

if ~isempty(length_adjusted_sessions)
    fprintf('Sessions with padded or size-adjusted label vectors: %s\n', mat2str(unique(length_adjusted_sessions(:)')));
end

if ~isempty(skipped_sessions)
    fprintf('Skipped entries: %s\n', mat2str(skipped_sessions(:)'));
end
end

function all_cells_path = resolve_all_cells_path()

repo_root = fullfile(fileparts(mfilename('fullpath')), '..', '..');
candidate_paths = { ...
    fullfile(repo_root, 'All_Cells_combined.mat'), ...
    fullfile(repo_root, 'Data', 'All_Cells_combined.mat'), ...
    fullfile(pwd, 'All_Cells_combined.mat')};

all_cells_path = first_existing_file(candidate_paths);
if isempty(all_cells_path)
    error('Could not find All_Cells_combined.mat. Pass all_cells_file or place it in the repository root/Data folder.')
end
end

function curated_session_root = resolve_curated_session_root()

repo_root = fullfile(fileparts(mfilename('fullpath')), '..', '..');
candidate_paths = { ...
    fullfile(repo_root, 'CCG_curated_connections', 'CuratedSessions'), ...
    fullfile(repo_root, 'Data', 'CCG_curated_connections', 'CuratedSessions'), ...
    fullfile(pwd, 'CCG_curated_connections', 'CuratedSessions')};

curated_session_root = first_existing_dir(candidate_paths);
if isempty(curated_session_root)
    curated_session_root = fullfile(repo_root, 'Data', 'CCG_curated_connections', 'CuratedSessions');
end
end

function path_out = first_existing_file(candidate_paths)

path_out = '';
for iPath = 1:numel(candidate_paths)
    if exist(candidate_paths{iPath}, 'file')
        path_out = candidate_paths{iPath};
        return
    end
end
end

function path_out = first_existing_dir(candidate_paths)

path_out = '';
for iPath = 1:numel(candidate_paths)
    if exist(candidate_paths{iPath}, 'dir')
        path_out = candidate_paths{iPath};
        return
    end
end
end

function unit_count = infer_all_cells_session_unit_count(session_entry)

unit_count = 0;

if isfield(session_entry, 'classific_firingRate') && ~isempty(session_entry.classific_firingRate)
    unit_count = numel(session_entry.classific_firingRate);
elseif isfield(session_entry, 'CCGbased_CellClassfication') && ~isempty(session_entry.CCGbased_CellClassfication)
    unit_count = numel(session_entry.CCGbased_CellClassfication);
end
end

function unit_count = infer_curated_session_unit_count(curated_session)

unit_count = 0;

if isfield(curated_session, 'n_units') && ~isempty(curated_session.n_units)
    unit_count = double(curated_session.n_units);
elseif isfield(curated_session, 'firing_rates') && ~isempty(curated_session.firing_rates)
    unit_count = numel(curated_session.firing_rates);
elseif isfield(curated_session, 'final_cell_classification') && ~isempty(curated_session.final_cell_classification)
    unit_count = numel(curated_session.final_cell_classification);
end
end
