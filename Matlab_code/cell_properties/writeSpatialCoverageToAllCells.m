function [All_Cells_combined, result] = writeSpatialCoverageToAllCells(allCellsInput, sessionIndex, behaviorName, rateMapInput, varargin)
% writeSpatialCoverageToAllCells calculates spatial coverage for one
% behavior and writes it into an All_Cells_combined session entry.
%
% allCellsInput can be either the All_Cells_combined struct array or a path
% to a MAT file containing All_Cells_combined. rateMapInput can be a rateMap
% struct, a MAT file containing rateMap, a cell array of maps, or a 3-D map
% stack with cells along the third dimension.
%
% Example:
%   [All_Cells_combined, result] = writeSpatialCoverageToAllCells( ...
%       fullfile(pwd, 'Data', 'All_Cells_combined.mat'), 11, 'of1', ...
%       fullfile(session_dir, 'processedData', 'PlaceMap', 'of1_Map.mat'), ...
%       'DoPlot', false);

opts = parse_writer_options(varargin{:});
if ~(ischar(behaviorName) || isstring(behaviorName))
    error('writeSpatialCoverageToAllCells:InvalidBehaviorName', ...
        'behaviorName must be a character vector or string scalar.');
end

[All_Cells_combined, all_cells_file, save_to_file, last_updated] = ...
    load_all_cells_input(allCellsInput, opts);

validateattributes(sessionIndex, {'numeric'}, {'scalar', 'integer', 'positive'});
if sessionIndex > numel(All_Cells_combined)
    error('writeSpatialCoverageToAllCells:SessionOutOfRange', ...
        'sessionIndex %d exceeds numel(All_Cells_combined) = %d.', ...
        sessionIndex, numel(All_Cells_combined));
end

[maps, cell_labels] = extract_rate_maps(rateMapInput);
coverage_values = NaN(numel(maps), 1);
plot_files = cell(numel(maps), 1);
field_name = resolve_output_field_name(behaviorName, opts.OutputFieldName);

for cellIdx = 1:numel(maps)
    plot_file = '';
    if opts.DoPlot
        plot_file = build_plot_file(opts, sessionIndex, behaviorName, cellIdx, cell_labels);
    end
    plot_files{cellIdx} = plot_file;
    close_plot = resolve_close_plot(opts, plot_file);

    coverage_values(cellIdx) = computeSpatialCoverage(maps{cellIdx}, ...
        'DoPlot', opts.DoPlot, ...
        'StopOnPlots', opts.StopOnPlots, ...
        'SaveFile', plot_file, ...
        'ClosePlot', close_plot, ...
        'FigureVisible', opts.FigureVisible, ...
        'Threshold', opts.Threshold, ...
        'Title', sprintf('%s cell %d', char(behaviorName), cellIdx));
end

All_Cells_combined(sessionIndex).(field_name) = coverage_values(:);

written_fields = {field_name};
if ~isempty(opts.MeanFieldName)
    source_fields = resolve_mean_source_fields(All_Cells_combined(sessionIndex), ...
        field_name, opts.MeanFieldName, opts.MeanSourceFieldNames);
    All_Cells_combined(sessionIndex).(opts.MeanFieldName) = ...
        column_mean_from_fields(All_Cells_combined(sessionIndex), source_fields);
    written_fields{end + 1} = opts.MeanFieldName;
end

if opts.UpdateLastUpdated
    last_updated = mark_spatial_coverage_fields_updated(last_updated, ...
        sessionIndex, written_fields, opts.UpdateScriptName);
end

if save_to_file
    save_all_cells_output(all_cells_file, All_Cells_combined, last_updated, opts);
end

result = struct( ...
    'sessionIndex', sessionIndex, ...
    'behaviorName', char(behaviorName), ...
    'outputFieldName', field_name, ...
    'coverageValues', coverage_values, ...
    'plotFiles', {plot_files}, ...
    'allCellsFile', all_cells_file);

if ~isempty(opts.MeanFieldName)
    result.meanFieldName = opts.MeanFieldName;
    result.meanSourceFieldNames = source_fields;
end
end


function opts = parse_writer_options(varargin)

opts = struct( ...
    'OutputFieldName', '', ...
    'MeanFieldName', '', ...
    'MeanSourceFieldNames', {{}}, ...
    'DoPlot', false, ...
    'StopOnPlots', false, ...
    'PlotSaveDir', '', ...
    'PlotFilePrefix', '', ...
    'ClosePlot', [], ...
    'FigureVisible', 'on', ...
    'Threshold', 0.75, ...
    'AllCellsFile', '', ...
    'SaveToFile', [], ...
    'MatFileVersion', '', ...
    'UpdateLastUpdated', true, ...
    'UpdateScriptName', mfilename);

if mod(numel(varargin), 2) ~= 0
    error('writeSpatialCoverageToAllCells:NameValueMissing', ...
        'Name-value options must be provided in pairs.');
end

for argIdx = 1:2:numel(varargin)
    option_name = normalize_option_name(varargin{argIdx});
    option_value = varargin{argIdx + 1};

    switch option_name
        case {'outputfieldname', 'outputfield'}
            opts.OutputFieldName = char(option_value);
        case {'meanfieldname', 'meanfield'}
            opts.MeanFieldName = char(option_value);
        case {'meansourcefieldnames', 'meansourcefields'}
            opts.MeanSourceFieldNames = cellstr(option_value);
        case {'doplot', 'plotting'}
            opts.DoPlot = to_logical(option_value);
        case {'stoponplots', 'pauseonplots'}
            opts.StopOnPlots = to_logical(option_value);
        case {'plotsavedir', 'plotdir'}
            opts.PlotSaveDir = char(option_value);
        case {'plotfileprefix', 'plotprefix'}
            opts.PlotFilePrefix = char(option_value);
        case {'closeplot', 'closeplots'}
            opts.ClosePlot = to_logical(option_value);
        case {'figurevisible', 'visible'}
            opts.FigureVisible = char(option_value);
        case 'threshold'
            opts.Threshold = double(option_value);
        case {'allcellsfile', 'savefile'}
            opts.AllCellsFile = char(option_value);
        case 'savetofile'
            opts.SaveToFile = to_logical(option_value);
        case {'matfileversion', 'saveversion'}
            opts.MatFileVersion = char(option_value);
        case {'updatelastupdated', 'trackupdates'}
            opts.UpdateLastUpdated = to_logical(option_value);
        case {'updatescriptname', 'updatescript'}
            opts.UpdateScriptName = char(option_value);
        otherwise
            error('writeSpatialCoverageToAllCells:UnknownOption', ...
                'Unknown option "%s".', char(varargin{argIdx}));
    end
end

if opts.StopOnPlots
    opts.DoPlot = true;
    opts.FigureVisible = 'on';
end

if opts.Threshold <= 0 || opts.Threshold > 1 || ~isfinite(opts.Threshold)
    error('writeSpatialCoverageToAllCells:InvalidThreshold', ...
        'Threshold must be finite and in the interval (0, 1].');
end
end


function [All_Cells_combined, all_cells_file, save_to_file, last_updated] = load_all_cells_input(allCellsInput, opts)

last_updated = struct();

if ischar(allCellsInput) || isstring(allCellsInput)
    all_cells_file = char(allCellsInput);
    loaded = load(all_cells_file, 'All_Cells_combined');
    if ~isfield(loaded, 'All_Cells_combined')
        error('writeSpatialCoverageToAllCells:MissingAllCells', ...
            'The file does not contain All_Cells_combined: %s', all_cells_file);
    end
    All_Cells_combined = loaded.All_Cells_combined;
    available_variables = who('-file', all_cells_file);
    if any(strcmp(available_variables, 'last_updated'))
        loaded_updates = load(all_cells_file, 'last_updated');
        last_updated = loaded_updates.last_updated;
    end
    default_save_to_file = true;
else
    All_Cells_combined = allCellsInput;
    all_cells_file = opts.AllCellsFile;
    if ~isempty(all_cells_file) && exist(all_cells_file, 'file') == 2
        available_variables = who('-file', all_cells_file);
        if any(strcmp(available_variables, 'last_updated'))
            loaded_updates = load(all_cells_file, 'last_updated');
            last_updated = loaded_updates.last_updated;
        end
    end
    default_save_to_file = ~isempty(all_cells_file);
end

if isempty(opts.SaveToFile)
    save_to_file = default_save_to_file;
else
    save_to_file = opts.SaveToFile;
end

if save_to_file && isempty(all_cells_file)
    error('writeSpatialCoverageToAllCells:MissingOutputFile', ...
        'SaveToFile is true, but no AllCellsFile was provided.');
end
end


function [maps, cell_labels] = extract_rate_maps(rateMapInput)

if ischar(rateMapInput) || isstring(rateMapInput)
    loaded = load(char(rateMapInput));
    if isfield(loaded, 'rateMap')
        rateMapInput = loaded.rateMap;
    else
        variable_names = fieldnames(loaded);
        if numel(variable_names) ~= 1
            error('writeSpatialCoverageToAllCells:MissingRateMap', ...
                'Could not identify the rateMap variable in %s.', char(rateMapInput));
        end
        rateMapInput = loaded.(variable_names{1});
    end
end

cell_labels = {};
if isstruct(rateMapInput)
    if ~isfield(rateMapInput, 'map')
        error('writeSpatialCoverageToAllCells:MissingMapField', ...
            'The rateMap struct must contain a map field.');
    end
    if isfield(rateMapInput, 'tFileNames')
        cell_labels = normalize_cell_labels(rateMapInput.tFileNames);
    end
    maps = normalize_maps(rateMapInput.map);
else
    maps = normalize_maps(rateMapInput);
end

if isempty(cell_labels)
    cell_labels = cell(numel(maps), 1);
end
cell_labels = cell_labels(:);
if numel(cell_labels) < numel(maps)
    cell_labels(end + 1:numel(maps), 1) = {''};
end
end


function cell_labels = normalize_cell_labels(cell_labels)

if isstring(cell_labels)
    cell_labels = cellstr(cell_labels);
elseif ischar(cell_labels)
    cell_labels = cellstr(cell_labels);
elseif ~iscell(cell_labels)
    cell_labels = {};
end
end


function maps = normalize_maps(mapInput)

if iscell(mapInput)
    maps = mapInput(:);
elseif isnumeric(mapInput)
    if ismatrix(mapInput)
        maps = {mapInput};
    elseif ndims(mapInput) == 3
        maps = cell(size(mapInput, 3), 1);
        for mapIdx = 1:size(mapInput, 3)
            maps{mapIdx} = mapInput(:, :, mapIdx);
        end
    else
        error('writeSpatialCoverageToAllCells:InvalidMapArray', ...
            'Numeric rate maps must be 2-D or 3-D.');
    end
else
    error('writeSpatialCoverageToAllCells:InvalidRateMapInput', ...
        'rateMapInput must provide a map cell array, 2-D map, or 3-D map stack.');
end

for mapIdx = 1:numel(maps)
    if ~isnumeric(maps{mapIdx}) || ndims(maps{mapIdx}) ~= 2
        error('writeSpatialCoverageToAllCells:InvalidMap', ...
            'Map %d is not a numeric 2-D matrix.', mapIdx);
    end
end
end


function field_name = resolve_output_field_name(behaviorName, outputFieldName)

if ~isempty(outputFieldName)
    field_name = outputFieldName;
    return
end

field_name = sprintf('spatial_coverage_%s', make_valid_field_part(behaviorName));
end


function plot_file = build_plot_file(opts, sessionIndex, behaviorName, cellIdx, cell_labels)

plot_file = '';
if isempty(opts.PlotSaveDir)
    return
end

if exist(opts.PlotSaveDir, 'dir') ~= 7
    mkdir(opts.PlotSaveDir);
end

if isempty(opts.PlotFilePrefix)
    plot_prefix = sprintf('SpatialCoverage_session%03d_%s', ...
        sessionIndex, make_valid_file_part(behaviorName));
else
    plot_prefix = opts.PlotFilePrefix;
end

cell_label = '';
if cellIdx <= numel(cell_labels)
    cell_label = cell_labels{cellIdx};
end

if isempty(cell_label)
    cell_part = sprintf('cell%03d', cellIdx);
else
    cell_part = make_valid_file_part(cell_label);
end

plot_file = fullfile(opts.PlotSaveDir, sprintf('%s_%s.tif', plot_prefix, cell_part));
end


function close_plot = resolve_close_plot(opts, plot_file)

if isempty(opts.ClosePlot)
    close_plot = opts.StopOnPlots || ~isempty(plot_file);
else
    close_plot = opts.ClosePlot;
end
end


function source_fields = resolve_mean_source_fields(session_data, output_field, mean_field, requested_fields)

if ~isempty(requested_fields)
    source_fields = requested_fields(:)';
    return
end

all_fields = fieldnames(session_data);
source_fields = {};
for fieldIdx = 1:numel(all_fields)
    current_field = all_fields{fieldIdx};
    is_spatial_coverage = strncmp(current_field, 'spatial_coverage_', numel('spatial_coverage_'));
    is_mean_field = strcmp(current_field, mean_field) || ~isempty(strfind(lower(current_field), 'mean'));
    if is_spatial_coverage && ~is_mean_field
        source_fields{end + 1} = current_field; %#ok<AGROW>
    end
end

if ~any(strcmp(source_fields, output_field))
    source_fields{end + 1} = output_field;
end
end


function mean_values = column_mean_from_fields(session_data, source_fields)

max_length = 0;
for fieldIdx = 1:numel(source_fields)
    values = get_session_column(session_data, source_fields{fieldIdx});
    max_length = max(max_length, numel(values));
end

coverage_matrix = NaN(max_length, numel(source_fields));
for fieldIdx = 1:numel(source_fields)
    values = get_session_column(session_data, source_fields{fieldIdx});
    coverage_matrix(1:numel(values), fieldIdx) = values(:);
end

valid_matrix = isfinite(coverage_matrix);
coverage_matrix(~valid_matrix) = 0;
valid_counts = sum(valid_matrix, 2);
mean_values = sum(coverage_matrix, 2) ./ valid_counts;
mean_values(valid_counts == 0) = NaN;
end


function values = get_session_column(session_data, field_name)

if ~isfield(session_data, field_name) || isempty(session_data.(field_name))
    values = [];
else
    values = double(session_data.(field_name)(:));
end
end


function last_updated = mark_spatial_coverage_fields_updated(last_updated, sessionIndex, field_names, updateScriptName)

if ~isstruct(last_updated)
    last_updated = struct();
end

timestamp = datestr(now, 'yyyy-mm-ddTHH:MM:SS');
for fieldIdx = 1:numel(field_names)
    last_updated(sessionIndex).(field_names{fieldIdx}) = struct( ...
        'update_timestamp', timestamp, ...
        'update_script', updateScriptName);
end
end


function save_all_cells_output(all_cells_file, All_Cells_combined, last_updated, opts)

[all_cells_dir, ~, ~] = fileparts(all_cells_file);
if ~isempty(all_cells_dir) && exist(all_cells_dir, 'dir') ~= 7
    mkdir(all_cells_dir);
end

if exist(all_cells_file, 'file') == 2
    save(all_cells_file, 'All_Cells_combined', '-append');
    if opts.UpdateLastUpdated
        save(all_cells_file, 'last_updated', '-append');
    end
elseif isempty(opts.MatFileVersion)
    if opts.UpdateLastUpdated
        save(all_cells_file, 'All_Cells_combined', 'last_updated');
    else
        save(all_cells_file, 'All_Cells_combined');
    end
else
    if opts.UpdateLastUpdated
        save(all_cells_file, 'All_Cells_combined', 'last_updated', opts.MatFileVersion);
    else
        save(all_cells_file, 'All_Cells_combined', opts.MatFileVersion);
    end
end
end


function field_part = make_valid_field_part(value)

if isstring(value)
    value = char(value);
end
field_part = regexprep(char(value), '[^A-Za-z0-9_]', '_');
field_part = regexprep(field_part, '_+', '_');
field_part = regexprep(field_part, '^_|_$', '');
if isempty(field_part)
    field_part = 'behavior';
end
if isempty(regexp(field_part(1), '[A-Za-z]', 'once'))
    field_part = ['x' field_part];
end
end


function file_part = make_valid_file_part(value)

if isstring(value)
    value = char(value);
end
file_part = regexprep(char(value), '[^A-Za-z0-9_.-]', '_');
file_part = regexprep(file_part, '_+', '_');
file_part = regexprep(file_part, '^_|_$', '');
if isempty(file_part)
    file_part = 'unnamed';
end
end


function option_name = normalize_option_name(option_name)

if isstring(option_name)
    option_name = char(option_name);
end
option_name = lower(strrep(char(option_name), '_', ''));
end


function value = to_logical(value)

if isempty(value)
    return
end

if islogical(value)
    value = any(value(:));
elseif isnumeric(value)
    value = any(value(:) ~= 0);
elseif isstring(value)
    value = char(value);
end

if ischar(value)
    value = any(strcmpi(value, {'true', 'on', 'yes', 'y', '1'}));
end
end
