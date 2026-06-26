function result = compute_place_cell_reuse_metrics(custom_settings)
% Compute place-cell reuse and Dice overlap between open-field contexts.

if nargin < 1 || isempty(custom_settings)
    custom_settings = struct();
end

all_cells_path = char(string(get_option_value(custom_settings, 'allCellsPath', '')));
if isempty(all_cells_path)
    all_cells_path = resolve_all_cells_path();
end

session_info_path = char(string(get_option_value(custom_settings, 'sessionInfoPath', '')));
if isempty(session_info_path)
    session_info_path = resolve_session_info_path();
end

session_indices = get_option_value(custom_settings, 'sessionIndices', []);
maze_prefix = char(string(get_option_value(custom_settings, 'mazePrefix', 'of')));
only_principal_cells = logical(get_option_value(custom_settings, 'onlyPrincipalCells', true));
include_optotagged_groups = logical(get_option_value(custom_settings, 'includeOptotaggedGroups', true));
save_updated_all_cells = logical(get_option_value(custom_settings, 'saveUpdatedAllCells', true));
save_session_files = logical(get_option_value(custom_settings, 'saveSessionFiles', true));
update_script_name = 'compute_place_cell_reuse_metrics.m';

loaded_all_cells = load(all_cells_path, 'All_Cells_combined');
if ~isfield(loaded_all_cells, 'All_Cells_combined')
    error('All_Cells_combined was not found in %s.', all_cells_path);
end
All_Cells_combined = loaded_all_cells.All_Cells_combined;
last_updated = load_last_updated(all_cells_path, numel(All_Cells_combined));

loaded_session_info = load(session_info_path, 'sessInfo');
if ~isfield(loaded_session_info, 'sessInfo')
    error('sessInfo was not found in %s.', session_info_path);
end
sessInfo = loaded_session_info.sessInfo;

if isempty(session_indices)
    session_indices = 1:min(numel(sessInfo), numel(All_Cells_combined));
end

pair_defs = open_field_pair_defs();
session_summaries = repmat(empty_reuse_summary(), numel(session_indices), 1);
summary_tables = {};

for idx = 1:numel(session_indices)
    session_index = session_indices(idx);
    session_summaries(idx).session_index = session_index;

    if session_index < 1 || session_index > numel(sessInfo) || session_index > numel(All_Cells_combined)
        session_summaries(idx).status = 'skipped_invalid_session_index';
        warning('Skipping invalid session index: %d', session_index);
        continue
    end

    session_data = All_Cells_combined(session_index);
    num_cells = infer_session_cell_count(session_data);
    if num_cells == 0
        session_summaries(idx).status = 'skipped_no_cell_fields';
        warning('Session %d: no cell-level fields available for place-cell reuse.', session_index);
        continue
    end

    map_availability = infer_open_field_map_availability(sessInfo(session_index), maze_prefix);
    [place_cell_code, of_available, code_source] = resolve_place_cell_context_code( ...
        session_data, num_cells, map_availability, only_principal_cells);

    if ~any(of_available)
        session_summaries(idx).status = 'skipped_no_open_field_contexts';
        warning('Session %d: no open-field place-cell contexts were available.', session_index);
        continue
    end

    [All_Cells_combined, written_fields, place_cell_reuse, session_summary_table] = ...
        write_session_reuse_fields(All_Cells_combined, session_index, place_cell_code, of_available, pair_defs, ...
        include_optotagged_groups);

    last_updated = mark_fields_updated(last_updated, session_index, written_fields, update_script_name);
    summary_tables{end + 1, 1} = session_summary_table; %#ok<AGROW>

    if save_session_files
        save_session_reuse_file(sessInfo(session_index), place_cell_code, of_available, place_cell_reuse, ...
            session_summary_table);
    end

    session_summaries(idx).status = 'processed';
    session_summaries(idx).cell_count = num_cells;
    session_summaries(idx).code_source = code_source;
    session_summaries(idx).of_available = of_available;
    session_summaries(idx).written_fields = written_fields;
end

reuse_summary_table = concatenate_tables(summary_tables);

if save_updated_all_cells
    save_updated_all_cells_file(all_cells_path, All_Cells_combined, last_updated);
end

result = struct( ...
    'allCellsPath', all_cells_path, ...
    'sessionInfoPath', session_info_path, ...
    'mazePrefix', maze_prefix, ...
    'sessions', session_summaries, ...
    'summaryTable', reuse_summary_table);

fprintf('Computed place-cell reuse metrics for %d sessions.\n', ...
    nnz(strcmp({session_summaries.status}, 'processed')));
end

function [All_Cells_combined, written_fields, place_cell_reuse, summary_table] = ...
        write_session_reuse_fields(All_Cells_combined, session_index, place_cell_code, of_available, pair_defs, ...
        include_optotagged_groups)

num_cells = size(place_cell_code, 1);
context_code = nan(num_cells, 3);
for context_idx = 1:3
    if of_available(context_idx)
        context_code(:, context_idx) = double(place_cell_code(:, context_idx) > 0);
    end
end

All_Cells_combined(session_index).place_cell_context_code = context_code;
All_Cells_combined(session_index).context_A_place_cell = context_code(:, 1);
All_Cells_combined(session_index).context_B_place_cell = context_code(:, 2);
All_Cells_combined(session_index).context_A_revisit_place_cell = context_code(:, 3);

written_fields = { ...
    'place_cell_context_code', ...
    'context_A_place_cell', ...
    'context_B_place_cell', ...
    'context_A_revisit_place_cell' ...
    };

place_cell_reuse = struct();
for pair_idx = 1:numel(pair_defs)
    pair_def = pair_defs(pair_idx);
    tested = false(num_cells, 1);
    reused = nan(num_cells, 1);

    if of_available(pair_def.contextAIndex) && of_available(pair_def.contextBIndex)
        tested(:) = true;
        reused(:) = double(place_cell_code(:, pair_def.contextAIndex) > 0 & ...
            place_cell_code(:, pair_def.contextBIndex) > 0);
    end

    reuse_field = sprintf('place_cell_reuse_%s', pair_def.fieldSuffix);
    tested_field = sprintf('place_cell_reuse_tested_%s', pair_def.fieldSuffix);

    All_Cells_combined(session_index).(reuse_field) = reused;
    All_Cells_combined(session_index).(tested_field) = double(tested);
    place_cell_reuse.(pair_def.fieldSuffix) = reused;
    place_cell_reuse.([pair_def.fieldSuffix '_tested']) = tested;

    written_fields{end + 1} = reuse_field; %#ok<AGROW>
    written_fields{end + 1} = tested_field; %#ok<AGROW>
end

summary_table = build_session_dice_summary( ...
    session_index, All_Cells_combined(session_index), place_cell_code, of_available, pair_defs, include_optotagged_groups);
All_Cells_combined(session_index).place_cell_reuse_dice_summary = summary_table;
written_fields{end + 1} = 'place_cell_reuse_dice_summary';
end

function summary_table = build_session_dice_summary( ...
        session_index, session_data, place_cell_code, of_available, pair_defs, include_optotagged_groups)

num_cells = size(place_cell_code, 1);
[animal_id, recording_day] = get_session_identifiers(session_data);
population_defs = build_population_defs(session_data, num_cells, include_optotagged_groups);
rows = cell(0, 13);

for population_idx = 1:size(population_defs, 1)
    population_label = population_defs{population_idx, 1};
    population_mask = population_defs{population_idx, 2};

    for pair_idx = 1:numel(pair_defs)
        pair_def = pair_defs(pair_idx);
        contexts_available = of_available(pair_def.contextAIndex) && of_available(pair_def.contextBIndex);
        valid_mask = false(num_cells, 1);
        if contexts_available
            valid_mask = population_mask(:);
        end

        positive_a = place_cell_code(:, pair_def.contextAIndex) > 0;
        positive_b = place_cell_code(:, pair_def.contextBIndex) > 0;

        context_a_n = sum(valid_mask & positive_a);
        context_b_n = sum(valid_mask & positive_b);
        both_n = sum(valid_mask & positive_a & positive_b);
        valid_cells_n = sum(valid_mask);

        if (context_a_n + context_b_n) > 0
            dice_value = (2 * both_n) / (context_a_n + context_b_n);
            dice_percent = 100 * dice_value;
        else
            dice_value = NaN;
            dice_percent = NaN;
        end

        rows(end + 1, :) = {session_index, animal_id, recording_day, population_label, ...
            pair_def.pairLabel, pair_def.contextAName, pair_def.contextBName, ...
            both_n, context_a_n, context_b_n, valid_cells_n, dice_value, dice_percent}; %#ok<AGROW>
    end
end

summary_table = cell2table(rows, 'VariableNames', ...
    {'SessionIndex', 'AnimalID', 'RecordingDay', 'Population', 'Pair', ...
    'ContextA', 'ContextB', 'Both_n', 'ContextA_n', 'ContextB_n', ...
    'ValidCells_n', 'Dice_value', 'Dice_percent'});
end

function population_defs = build_population_defs(session_data, num_cells, include_optotagged_groups)

population_defs = {'all', true(num_cells, 1)};
if ~include_optotagged_groups
    return
end

optotagged = normalize_numeric_cell_vector(get_first_available_field(session_data, {'optotagged'}), num_cells, NaN);
if isempty(optotagged) || all(isnan(optotagged))
    return
end

population_defs(end + 1, :) = {'cfos-', optotagged <= 0}; %#ok<AGROW>
population_defs(end + 1, :) = {'cfos+', optotagged == 1}; %#ok<AGROW>
end

function [place_cell_code, of_available, code_source] = resolve_place_cell_context_code( ...
        session_data, num_cells, map_availability, only_principal_cells)

place_cell_code = zeros(num_cells, 3);
of_available = logical(map_availability(:)');
code_source = 'empty';

raw_code = get_first_available_field(session_data, {'place_cell_OF_combinations_code', 'place_cell_context_code'});
if ~isempty(raw_code)
    normalized_code = normalize_context_code(raw_code, num_cells);
    if ~isempty(normalized_code)
        place_cell_code = double(normalized_code > 0);
        code_source = 'place_cell_OF_combinations_code';
        of_available = of_available | any(isfinite(normalized_code), 1);
    end
end

principal_mask = infer_principal_cell_mask(session_data, num_cells, only_principal_cells);
field_sources = { ...
    {'PF_fieldNumbers_of1', 'of1_place_field_numbers'}, ...
    {'PF_fieldNumbers_of2', 'of2_place_field_numbers'}, ...
    {'PF_fieldNumbers_of3', 'of3_place_field_numbers'} ...
    };

used_place_field_numbers = false;
for context_idx = 1:3
    place_field_counts = normalize_numeric_cell_vector( ...
        get_first_available_field(session_data, field_sources{context_idx}), num_cells, NaN);
    if isempty(place_field_counts)
        continue
    end
    of_available(context_idx) = true;
    place_cell_code(:, context_idx) = double(principal_mask & place_field_counts > 0);
    used_place_field_numbers = true;
end

if used_place_field_numbers
    if strcmp(code_source, 'empty')
        code_source = 'place_field_numbers';
    else
        code_source = [code_source '+place_field_numbers'];
    end
end

for context_idx = 1:3
    if ~of_available(context_idx)
        place_cell_code(:, context_idx) = 0;
    end
end
end

function principal_mask = infer_principal_cell_mask(session_data, num_cells, only_principal_cells)

principal_mask = true(num_cells, 1);
if ~only_principal_cells
    return
end

if isfield(session_data, 'final_classification_numeric') && ~isempty(session_data.final_classification_numeric)
    classification = normalize_numeric_cell_vector(session_data.final_classification_numeric, num_cells, NaN);
    principal_mask = ismember(classification, [1 2]);
elseif isfield(session_data, 'GMM_based_classification_days') && ~isempty(session_data.GMM_based_classification_days)
    classification = normalize_numeric_cell_vector(session_data.GMM_based_classification_days, num_cells, NaN);
    principal_mask = classification == 1;
end
end

function context_code = normalize_context_code(raw_code, num_cells)

context_code = [];
raw_code = double(raw_code);

if isempty(raw_code)
    return
end

if size(raw_code, 1) == num_cells && size(raw_code, 2) >= 3
    context_code = raw_code(:, 1:3);
elseif size(raw_code, 2) == num_cells && size(raw_code, 1) >= 3
    context_code = raw_code(1:3, :)';
end
end

function num_cells = infer_session_cell_count(session_data)

num_cells = 0;
candidate_fields = { ...
    'final_classification_numeric', ...
    'GMM_based_classification_days', ...
    'optotagged', ...
    'of1_place_field_numbers', ...
    'of2_place_field_numbers', ...
    'of3_place_field_numbers', ...
    'PF_fieldNumbers_of1', ...
    'PF_fieldNumbers_of2', ...
    'PF_fieldNumbers_of3', ...
    'context_A_place_cell', ...
    'context_B_place_cell', ...
    'context_A_revisit_place_cell' ...
    };

for field_idx = 1:numel(candidate_fields)
    values = get_first_available_field(session_data, candidate_fields(field_idx));
    if ~isempty(values)
        num_cells = max(num_cells, numel(values));
    end
end

raw_code = get_first_available_field(session_data, {'place_cell_OF_combinations_code', 'place_cell_context_code'});
if ~isempty(raw_code)
    raw_size = size(raw_code);
    if numel(raw_size) >= 2
        if raw_size(2) == 3
            num_cells = max(num_cells, raw_size(1));
        elseif raw_size(1) == 3
            num_cells = max(num_cells, raw_size(2));
        end
    end
end
end

function availability = infer_open_field_map_availability(sess_entry, maze_prefix)

availability = false(1, 3);
if ~isfield(sess_entry, 'mainDir') || isempty(sess_entry.mainDir)
    return
end

place_map_dir = fullfile(sess_entry.mainDir, 'processedData', 'PlaceMap');
for context_idx = 1:3
    map_path = fullfile(place_map_dir, sprintf('%s%d_Map.mat', maze_prefix, context_idx));
    availability(context_idx) = exist(map_path, 'file') == 2;
end
end

function pair_defs = open_field_pair_defs()

pair_defs = struct( ...
    'fieldSuffix', {}, ...
    'pairLabel', {}, ...
    'contextAName', {}, ...
    'contextBName', {}, ...
    'contextAIndex', {}, ...
    'contextBIndex', {});

pair_defs(1).fieldSuffix = 'context_A_context_B';
pair_defs(1).pairLabel = 'Context A-Context B';
pair_defs(1).contextAName = 'Context A';
pair_defs(1).contextBName = 'Context B';
pair_defs(1).contextAIndex = 1;
pair_defs(1).contextBIndex = 2;

pair_defs(2).fieldSuffix = 'context_B_context_A_revisit';
pair_defs(2).pairLabel = 'Context B-Context A Revisit';
pair_defs(2).contextAName = 'Context B';
pair_defs(2).contextBName = 'Context A Revisit';
pair_defs(2).contextAIndex = 2;
pair_defs(2).contextBIndex = 3;

pair_defs(3).fieldSuffix = 'context_A_context_A_revisit';
pair_defs(3).pairLabel = 'Context A-Context A Revisit';
pair_defs(3).contextAName = 'Context A';
pair_defs(3).contextBName = 'Context A Revisit';
pair_defs(3).contextAIndex = 1;
pair_defs(3).contextBIndex = 3;
end

function save_session_reuse_file(sess_entry, place_cell_code, of_available, place_cell_reuse, summary_table)

if ~isfield(sess_entry, 'mainDir') || isempty(sess_entry.mainDir)
    return
end

output_dir = fullfile(sess_entry.mainDir, 'processedData', 'PlaceCellReuse');
if exist(output_dir, 'dir') ~= 7
    mkdir(output_dir);
end

placeCellContextCode = place_cell_code; %#ok<NASGU>
openFieldAvailable = of_available; %#ok<NASGU>
placeCellReuse = place_cell_reuse; %#ok<NASGU>
placeCellReuseSummary = summary_table; %#ok<NASGU>
save(fullfile(output_dir, 'place_cell_reuse_metrics.mat'), ...
    'placeCellContextCode', 'openFieldAvailable', 'placeCellReuse', 'placeCellReuseSummary');

try
    writetable(summary_table, fullfile(output_dir, 'place_cell_reuse_dice_summary.csv'));
catch ME
    warning('Could not write place-cell reuse CSV for %s: %s', output_dir, ME.message);
end
end

function [animal_id, recording_day] = get_session_identifiers(session_data)

animal_id = 'missing';
recording_day = NaN;

if isfield(session_data, 'animal') && ~isempty(session_data.animal)
    animal_id = normalize_label(session_data.animal);
end
if isfield(session_data, 'rec_day') && ~isempty(session_data.rec_day)
    recording_day = normalize_numeric_scalar(session_data.rec_day);
end
end

function label = normalize_label(value)

if iscell(value)
    if isempty(value)
        label = 'missing';
        return
    end
    value = value{1};
end

if isempty(value)
    label = 'missing';
elseif isnumeric(value) || islogical(value)
    label = num2str(value(1));
elseif isstring(value)
    label = char(value(1));
elseif ischar(value)
    if size(value, 1) > 1
        label = strtrim(value(1, :));
    else
        label = strtrim(value);
    end
else
    label = class(value);
end
end

function value = normalize_numeric_scalar(values)

values = double(values(:));
if isempty(values)
    value = NaN;
else
    value = values(1);
end
end

function values = normalize_numeric_cell_vector(values, num_cells, default_value)

if isempty(values)
    values = [];
    return
end

if nargin < 3 || isempty(default_value)
    default_value = NaN;
end

raw_values = double(values(:));
values = repmat(default_value, num_cells, 1);
n = min(num_cells, numel(raw_values));
values(1:n) = raw_values(1:n);
end

function values = get_first_available_field(session_data, field_names)

values = [];
if ischar(field_names) || isstring(field_names)
    field_names = cellstr(string(field_names));
end

for field_idx = 1:numel(field_names)
    field_name = field_names{field_idx};
    if isfield(session_data, field_name) && ~isempty(session_data.(field_name))
        values = session_data.(field_name);
        return
    end
end
end

function table_out = concatenate_tables(table_cells)

if isempty(table_cells)
    table_out = cell2table(cell(0, 13), 'VariableNames', ...
        {'SessionIndex', 'AnimalID', 'RecordingDay', 'Population', 'Pair', ...
        'ContextA', 'ContextB', 'Both_n', 'ContextA_n', 'ContextB_n', ...
        'ValidCells_n', 'Dice_value', 'Dice_percent'});
    return
end

table_out = table_cells{1};
for table_idx = 2:numel(table_cells)
    table_out = [table_out; table_cells{table_idx}]; %#ok<AGROW>
end
end

function summary = empty_reuse_summary()

summary = struct( ...
    'session_index', NaN, ...
    'status', '', ...
    'cell_count', 0, ...
    'code_source', '', ...
    'of_available', false(1, 3), ...
    'written_fields', {{}} ...
    );
end

function all_cells_path = resolve_all_cells_path()

repo_root = fullfile(fileparts(mfilename('fullpath')), '..', '..');
candidate_paths = { ...
    fullfile(repo_root, 'All_Cells_combined.mat'), ...
    fullfile(repo_root, 'Data', 'All_Cells_combined.mat'), ...
    fullfile(pwd, 'All_Cells_combined.mat')};

all_cells_path = first_existing_file(candidate_paths);
if isempty(all_cells_path)
    error('Could not find All_Cells_combined.mat. Pass custom_settings.allCellsPath or place it in the repository root/Data folder.');
end
end

function session_info_path = resolve_session_info_path()

repo_root = fullfile(fileparts(mfilename('fullpath')), '..', '..');
candidate_paths = { ...
    fullfile(repo_root, 'sessionInfo.mat'), ...
    fullfile(repo_root, 'Data', 'sessionInfo.mat'), ...
    fullfile(repo_root, 'Analysis_scripts', 'DataOrganization', 'sessionInfo.mat'), ...
    fullfile(pwd, 'sessionInfo.mat')};

session_info_path = first_existing_file(candidate_paths);
if isempty(session_info_path)
    error('Could not find sessionInfo.mat. Pass custom_settings.sessionInfoPath or place it in the repository root/Data folder.');
end
end

function path_out = first_existing_file(candidate_paths)

path_out = '';
for path_idx = 1:numel(candidate_paths)
    if exist(candidate_paths{path_idx}, 'file') == 2
        path_out = candidate_paths{path_idx};
        return
    end
end
end

function last_updated = load_last_updated(all_cells_path, session_count)

last_updated = struct();
if exist(all_cells_path, 'file') ~= 2
    return
end

available_variables = who('-file', all_cells_path);
if any(strcmp(available_variables, 'last_updated'))
    loaded = load(all_cells_path, 'last_updated');
    if isstruct(loaded.last_updated)
        last_updated = loaded.last_updated;
    end
end

last_updated = ensure_last_updated_session_capacity(last_updated, session_count);
end

function save_updated_all_cells_file(all_cells_path, All_Cells_combined, last_updated)

if exist(all_cells_path, 'file') == 2
    save(all_cells_path, 'All_Cells_combined', 'last_updated', '-append');
else
    save(all_cells_path, 'All_Cells_combined', 'last_updated', '-v7.3');
end
end

function last_updated = mark_fields_updated(last_updated, session_index, fields, update_script_name)

last_updated = ensure_last_updated_session_capacity(last_updated, session_index);
timestamp = datestr(now, 'yyyy-mm-ddTHH:MM:SS');
for field_idx = 1:numel(fields)
    field_name = fields{field_idx};
    if ~isfield(last_updated.sessions(session_index).fields, field_name)
        last_updated.sessions(session_index).fields.(field_name) = struct();
    end
    last_updated.sessions(session_index).fields.(field_name).timestamp = timestamp;
    last_updated.sessions(session_index).fields.(field_name).script = update_script_name;
end
end

function last_updated = ensure_last_updated_session_capacity(last_updated, session_count)

if ~isstruct(last_updated) || ~isfield(last_updated, 'sessions')
    last_updated = struct('sessions', struct('fields', {}));
end

if numel(last_updated.sessions) < session_count
    current_count = numel(last_updated.sessions);
    for session_idx = current_count + 1:session_count
        last_updated.sessions(session_idx).fields = struct();
    end
end

for session_idx = 1:session_count
    if ~isfield(last_updated.sessions(session_idx), 'fields') || isempty(last_updated.sessions(session_idx).fields)
        last_updated.sessions(session_idx).fields = struct();
    end
end
end

function value = get_option_value(settings_struct, field_name, default_value)

value = default_value;
if nargin < 1 || isempty(settings_struct) || ~isstruct(settings_struct)
    return
end

if isfield(settings_struct, field_name) && ~isempty(settings_struct.(field_name))
    value = settings_struct.(field_name);
end
end
