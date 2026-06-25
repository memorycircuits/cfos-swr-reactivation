function [spatial_coverage, details] = computeSpatialCoverage(rate_map, varargin)
% M Haberl, 05/04/2025
% Last edited: 2026-05-13
%
% 2D implementation to compute the spatial coverage metric of neurons.
% The metric is based on the higher coverage of SOM+ and PV+ interneurons
% compared to pyramidal neurons, see Royer et al. 2012, Supplementary
% Figure 7B. Royer calculated coverage from linear maps; this function
% adapts that approach to 2D rate maps. PYR cells are typically below 0.6,
% whereas INT cells are typically above 0.6.
%
% computeSpatialCoverage calculates the fraction of visited bins required
% to explain a fixed fraction of total map activity.
%
% The default threshold is 0.75, matching the previous implementation.
% Plotting is off by default. The previous plotting struct input is still
% supported:
%   plotting.do_plot = 1;
%   plotting.savefile_coverage = 'coverage_plot.tif';
%   spatial_coverage = computeSpatialCoverage(rate_map, plotting);

opts = parse_spatial_coverage_options(varargin{:});

rates = double(rate_map(:));
rates = rates(isfinite(rates));

sorted_rates = sort(rates, 'descend');
cum_sum = cumsum(sorted_rates);
valid_bin_count = numel(rates);
bins_to_threshold = NaN;
norm_cum_sum = zeros(size(cum_sum));

if valid_bin_count == 0
    spatial_coverage = NaN;
else
    total_activity = max(cum_sum);
    if total_activity <= 0 || ~isfinite(total_activity)
        bins_to_threshold = 0;
        spatial_coverage = 0;
    else
        norm_cum_sum = cum_sum ./ total_activity;
        bins_to_threshold = find(norm_cum_sum >= opts.Threshold, 1, 'first');
        if isempty(bins_to_threshold)
            bins_to_threshold = 0;
        end
        spatial_coverage = bins_to_threshold / valid_bin_count;
    end
end

details = struct( ...
    'threshold', opts.Threshold, ...
    'valid_bin_count', valid_bin_count, ...
    'bins_to_threshold', bins_to_threshold, ...
    'sorted_rates', sorted_rates, ...
    'norm_cum_sum', norm_cum_sum);

if opts.DoPlot
    plot_spatial_coverage(rate_map, sorted_rates, norm_cum_sum, ...
        bins_to_threshold, spatial_coverage, opts);
end
end


function opts = parse_spatial_coverage_options(varargin)

opts = struct( ...
    'DoPlot', false, ...
    'StopOnPlots', false, ...
    'SaveFile', '', ...
    'ClosePlot', [], ...
    'FigureVisible', 'on', ...
    'Threshold', 0.75, ...
    'Title', '');

arg_idx = 1;
if ~isempty(varargin) && isstruct(varargin{1})
    opts = apply_spatial_coverage_option_struct(opts, varargin{1});
    arg_idx = 2;
end

while arg_idx <= numel(varargin)
    if arg_idx == numel(varargin)
        error('computeSpatialCoverage:NameValueMissing', ...
            'Name-value options must be provided in pairs.');
    end

    option_name = normalize_option_name(varargin{arg_idx});
    option_value = varargin{arg_idx + 1};

    switch option_name
        case {'doplot', 'plotting'}
            opts.DoPlot = to_logical(option_value);
        case {'stoponplots', 'pauseonplots'}
            opts.StopOnPlots = to_logical(option_value);
        case {'savefile', 'savefilecoverage', 'plotfile'}
            opts.SaveFile = char(option_value);
        case {'closeplot', 'closeplots'}
            opts.ClosePlot = to_logical(option_value);
        case {'figurevisible', 'visible'}
            opts.FigureVisible = char(option_value);
        case 'threshold'
            opts.Threshold = double(option_value);
        case {'title', 'plottitle'}
            opts.Title = char(option_value);
        otherwise
            error('computeSpatialCoverage:UnknownOption', ...
                'Unknown option "%s".', char(varargin{arg_idx}));
    end

    arg_idx = arg_idx + 2;
end

if opts.Threshold <= 0 || opts.Threshold > 1 || ~isfinite(opts.Threshold)
    error('computeSpatialCoverage:InvalidThreshold', ...
        'Threshold must be finite and in the interval (0, 1].');
end

if opts.StopOnPlots
    opts.DoPlot = true;
    opts.FigureVisible = 'on';
end

if isempty(opts.ClosePlot)
    opts.ClosePlot = opts.StopOnPlots || ~isempty(opts.SaveFile);
end
end


function opts = apply_spatial_coverage_option_struct(opts, option_struct)

field_names = fieldnames(option_struct);
for iField = 1:numel(field_names)
    field_name = field_names{iField};
    option_name = normalize_option_name(field_name);
    option_value = option_struct.(field_name);

    switch option_name
        case {'doplot', 'plotting'}
            opts.DoPlot = to_logical(option_value);
        case {'stoponplots', 'pauseonplots', 'stoponplot'}
            opts.StopOnPlots = to_logical(option_value);
        case {'savefilecoverage', 'savefile', 'plotfile'}
            opts.SaveFile = char(option_value);
        case {'closeplot', 'closeplots'}
            opts.ClosePlot = to_logical(option_value);
        case {'figurevisible', 'visible'}
            opts.FigureVisible = char(option_value);
        case 'threshold'
            opts.Threshold = double(option_value);
        case {'title', 'plottitle'}
            opts.Title = char(option_value);
    end
end
end


function plot_spatial_coverage(rate_map, sorted_rates, norm_cum_sum, ...
        bins_to_threshold, spatial_coverage, opts)

fig = figure('Visible', opts.FigureVisible);

if ~isempty(opts.Title)
    set(fig, 'Name', opts.Title, 'NumberTitle', 'off');
end

subplot(2, 2, 1);
imagesc(rate_map);
axis image
axis off
axis equal
set(gca, 'YDir', 'normal');
colorbar;
title('Original Rate Map');

subplot(2, 2, 2);
bar(sorted_rates, 'k');
title('Bins Sorted by Rate');
ylabel('Hz');

subplot(2, 2, [3 4]);
plot(norm_cum_sum * 100, 'b', 'LineWidth', 2);
hold on
yline(opts.Threshold * 100, 'k--');
if isfinite(bins_to_threshold)
    xline(bins_to_threshold, 'k--');
end
xlabel('Bins (sorted)');
ylabel('%');
title('Normalized Cumulative Sum');
legend('Cumulative Sum', 'Threshold', 'Location', 'southeast');

if isfinite(bins_to_threshold)
    text(bins_to_threshold + 5, 50, sprintf('Coverage = %.2f', spatial_coverage));
else
    text(1, 50, sprintf('Coverage = %.2f', spatial_coverage));
end

if ~isempty(opts.SaveFile)
    [save_dir, ~, ~] = fileparts(opts.SaveFile);
    if ~isempty(save_dir) && exist(save_dir, 'dir') ~= 7
        mkdir(save_dir);
    end
    saveas(fig, opts.SaveFile);
end

if opts.StopOnPlots
    add_continue_button(fig);
    uiwait(fig);
    if ishandle(fig)
        set(fig, 'CloseRequestFcn', 'closereq');
    end
end

if opts.ClosePlot && ishandle(fig)
    close(fig);
end
end


function add_continue_button(fig)

set(fig, 'CloseRequestFcn', @(src, ~) uiresume(src));
uicontrol(fig, ...
    'Style', 'pushbutton', ...
    'String', 'Continue', ...
    'Units', 'normalized', ...
    'Position', [0.82 0.02 0.14 0.05], ...
    'Callback', @(~, ~) uiresume(fig));
end


function option_name = normalize_option_name(option_name)

if isstring(option_name)
    option_name = char(option_name);
end
option_name = lower(strrep(option_name, '_', ''));
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

