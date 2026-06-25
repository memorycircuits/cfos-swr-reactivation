function batch_export_summary = review_CCG_connection_labels(review_input, process_mode, export_options)

batch_export_summary = [];

if nargin < 3 || isempty(export_options)
    export_options = struct();
end

path_options = normalize_path_options(export_options);
config = load_classification_config(get_option_value(path_options, 'configPath', ''));
add_repo_matlab_code_path();
add_configured_dependency_paths(config, path_options);

review_root = char(string(get_option_value(path_options, 'reviewRoot', '')));
if isempty(review_root)
    review_root = default_review_root();
end

all_cells_file = char(string(get_option_value(path_options, 'allCellsPath', '')));
if isempty(all_cells_file)
    all_cells_file = resolve_all_cells_path();
end

curated_root = char(string(get_option_value(path_options, 'curatedRoot', '')));
if isempty(curated_root)
    curated_root = default_curated_root();
end

curated_session_root = char(string(get_option_value(path_options, 'curatedSessionRoot', '')));
if isempty(curated_session_root)
    curated_session_root = fullfile(curated_root, 'CuratedSessions');
end

summary_root = char(string(get_option_value(path_options, 'summaryRoot', '')));
if isempty(summary_root)
    summary_root = fullfile(curated_root, 'Summaries');
end

cell_pairs_root = char(string(get_option_value(path_options, 'cellPairsRoot', '')));
if isempty(cell_pairs_root)
    cell_pairs_root = default_cell_pairs_root();
end

accepted_plot_root = char(string(get_option_value(path_options, 'acceptedPlotRoot', '')));
if isempty(accepted_plot_root)
    accepted_plot_root = fullfile(cell_pairs_root, 'AcceptedCurated');
end

ensure_dir(curated_root);
ensure_dir(curated_session_root);
ensure_dir(summary_root);
ensure_dir(cell_pairs_root);
ensure_dir(accepted_plot_root);

if nargin < 1 || isempty(review_input)
    [review_files, process_mode] = select_review_targets(review_root);
else
    if nargin < 2 || isempty(process_mode)
        process_mode = 'only_new';
    end
    review_files = resolve_review_files(review_input, review_root);
end

if isempty(review_files)
    return
end

if strcmp(process_mode, 'export_plots_only')
    batch_export_summary = export_accepted_connection_plots_from_review_files( ...
        review_files, accepted_plot_root, export_options);
    return
end

[review_files, session_cache, skipped_mode_count, skipped_empty_count, skipped_empty_summary] = ...
    filter_review_files(review_files, process_mode);
if isempty(review_files)
    msgbox(sprintf(['No review files matched the selected mode.\n' ...
        'Skipped by mode: %d\nSkipped with no candidate connections: %d'], ...
        skipped_mode_count, skipped_empty_count), ...
        'CCG Review');
    return
end

animal_groups = build_animal_groups(session_cache);
animal_labels = build_animal_labels(session_cache, animal_groups);
current_animal_idx = 1;
active_file_indices = animal_groups{current_animal_idx};

current_file_idx = 1;
current_idx = 1;
current_review_file = '';
review_session = struct();
session_dirty = false(numel(review_files), 1);
presyn_fig = [];

fig = figure( ...
    'Name', 'CCG Connection Review', ...
    'Units', 'normalized', ...
    'Position', [0.03 0.05 0.94 0.90], ...
    'Color', 'w', ...
    'NumberTitle', 'off', ...
    'CloseRequestFcn', @close_reviewer ...
    );

txt_mode = uicontrol( ...
    'Parent', fig, ...
    'Style', 'text', ...
    'Units', 'normalized', ...
    'Position', [0.04 0.95 0.92 0.03], ...
    'BackgroundColor', 'w', ...
    'HorizontalAlignment', 'left', ...
    'FontSize', 10, ...
    'String', '' ...
    );

popup_session = uicontrol( ...
    'Parent', fig, ...
    'Style', 'popupmenu', ...
    'Units', 'normalized', ...
    'Position', [0.04 0.91 0.44 0.035], ...
    'String', {'Loading...'}, ...
    'Callback', @jump_session ...
    );

popup_pair = uicontrol( ...
    'Parent', fig, ...
    'Style', 'popupmenu', ...
    'Units', 'normalized', ...
    'Position', [0.52 0.91 0.44 0.035], ...
    'String', {'Loading...'}, ...
    'Callback', @jump_pair ...
    );

txt_info = uicontrol( ...
    'Parent', fig, ...
    'Style', 'text', ...
    'Units', 'normalized', ...
    'Position', [0.04 0.81 0.92 0.08], ...
    'BackgroundColor', 'w', ...
    'HorizontalAlignment', 'left', ...
    'FontSize', 11, ...
    'Max', 3, ...
    'String', '' ...
    );

ax_full = axes('Parent', fig, 'Units', 'normalized', 'Position', [0.04 0.35 0.28 0.40]);
ax_zoom = axes('Parent', fig, 'Units', 'normalized', 'Position', [0.36 0.35 0.28 0.40]);
ax_raw = axes('Parent', fig, 'Units', 'normalized', 'Position', [0.68 0.35 0.28 0.40]);

btn_prev_session = uicontrol('Parent', fig, 'Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', [0.04 0.26 0.11 0.05], 'String', 'Prev Session', ...
    'Callback', @go_previous_session);
btn_next_session = uicontrol('Parent', fig, 'Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', [0.16 0.26 0.11 0.05], 'String', 'Next Session', ...
    'Callback', @go_next_session);
btn_next_animal = uicontrol('Parent', fig, 'Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', [0.28 0.11 0.16 0.05], 'String', 'Proceed to Next Animal', ...
    'Callback', @go_next_animal);
btn_prev_pair = uicontrol('Parent', fig, 'Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', [0.31 0.26 0.11 0.05], 'String', 'Previous Pair', ...
    'Callback', @go_previous_pair);
btn_next_pair = uicontrol('Parent', fig, 'Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', [0.43 0.26 0.11 0.05], 'String', 'Next Pair', ...
    'Callback', @go_next_pair);

btn_accept = uicontrol('Parent', fig, 'Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', [0.04 0.18 0.12 0.05], 'String', 'Accept Auto', ...
    'Callback', @accept_auto);
btn_discard = uicontrol('Parent', fig, 'Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', [0.17 0.18 0.10 0.05], 'String', 'Discard', ...
    'Callback', @discard_pair);
btn_switch = uicontrol('Parent', fig, 'Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', [0.28 0.18 0.10 0.05], 'String', 'Switch Type', ...
    'Callback', @switch_type);
btn_set_exc = uicontrol('Parent', fig, 'Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', [0.39 0.18 0.10 0.05], 'String', 'Set Excit.', ...
    'Callback', @set_excitatory);
btn_set_inh = uicontrol('Parent', fig, 'Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', [0.50 0.18 0.10 0.05], 'String', 'Set Inhib.', ...
    'Callback', @set_inhibitory);
btn_save = uicontrol('Parent', fig, 'Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', [0.64 0.18 0.11 0.05], 'String', 'Save Session', ...
    'Callback', @save_review);
btn_write = uicontrol('Parent', fig, 'Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', [0.76 0.18 0.17 0.05], 'String', 'Write Session Labels', ...
    'Callback', @write_session_labels);
btn_global_summary = uicontrol('Parent', fig, 'Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', [0.76 0.11 0.17 0.05], 'String', 'Write Global Summary', ...
    'Callback', @write_global_summary_callback);
btn_examine_pre = uicontrol('Parent', fig, 'Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', [0.04 0.11 0.22 0.05], 'String', 'Examine Presynaptic Cell', ...
    'Callback', @examine_presynaptic_cell);

txt_status = uicontrol( ...
    'Parent', fig, ...
    'Style', 'text', ...
    'Units', 'normalized', ...
    'Position', [0.04 0.02 0.92 0.07], ...
    'BackgroundColor', 'w', ...
    'HorizontalAlignment', 'left', ...
    'FontSize', 10, ...
    'Max', 2, ...
    'String', '' ...
    );

set(txt_mode, 'String', build_mode_string(process_mode, numel(review_files), ...
    skipped_mode_count, skipped_empty_count, animal_labels{current_animal_idx}, ...
    current_animal_idx, numel(animal_groups), numel(active_file_indices)));

load_session(active_file_indices(1), true);

    function go_previous_session(~, ~)
        local_idx = current_active_session_position();
        if local_idx > 1
            load_session(active_file_indices(local_idx - 1), true);
        end
    end

    function go_next_session(~, ~)
        local_idx = current_active_session_position();
        if local_idx < numel(active_file_indices)
            load_session(active_file_indices(local_idx + 1), true);
        end
    end

    function jump_session(src, ~)
        local_idx = get(src, 'Value');
        target_idx = active_file_indices(local_idx);
        if target_idx ~= current_file_idx
            load_session(target_idx, true);
        end
    end

    function go_next_animal(~, ~)
        if current_animal_idx >= numel(animal_groups)
            return
        end
        switch_to_animal(current_animal_idx + 1);
    end

    function go_previous_pair(~, ~)
        if isempty(review_session.review_pairs)
            return
        end
        current_idx = max(1, current_idx - 1);
        render_pair();
    end

    function go_next_pair(~, ~)
        if isempty(review_session.review_pairs)
            return
        end
        current_idx = min(numel(review_session.review_pairs), current_idx + 1);
        render_pair();
    end

    function jump_pair(src, ~)
        if isempty(review_session.review_pairs)
            return
        end
        current_idx = get(src, 'Value');
        render_pair();
    end

    function accept_auto(~, ~)
        if isempty(review_session.review_pairs)
            return
        end
        apply_decision(review_session.review_pairs(current_idx).initial_decision);
    end

    function discard_pair(~, ~)
        if isempty(review_session.review_pairs)
            return
        end
        apply_decision(0);
    end

    function switch_type(~, ~)
        if isempty(review_session.review_pairs)
            return
        end

        decision = review_session.review_pairs(current_idx).current_decision;
        if decision == 1
            decision = 2;
        elseif decision == 2
            decision = 1;
        elseif review_session.review_pairs(current_idx).initial_decision == 1
            decision = 2;
        else
            decision = 1;
        end

        apply_decision(decision);
    end

    function set_excitatory(~, ~)
        if isempty(review_session.review_pairs)
            return
        end
        apply_decision(1);
    end

    function set_inhibitory(~, ~)
        if isempty(review_session.review_pairs)
            return
        end
        apply_decision(2);
    end

    function apply_decision(decision_value)
        review_session.review_pairs(current_idx).current_decision = decision_value;
        review_session.review_pairs(current_idx).verified = true;
        review_session.review_pairs(current_idx).last_reviewed = char(datetime('now'));
        review_session.writeback_complete = false;
        review_session.last_writeback = '';
        review_session.curated_file = '';
        review_session.accepted_plot_folder = '';
        review_session.accepted_plot_count = 0;
        review_session = synchronize_review_session_status(review_session);
        sync_current_session_cache();
        session_dirty(current_file_idx) = true;
        refresh_selectors();
        render_pair();
    end

    function save_review(~, ~)
        save_current_session(true);
        set(txt_status, 'String', sprintf('Saved review session: %s', current_review_file));
    end

    function write_session_labels(~, ~)
        set(txt_status, 'String', sprintf('Started writing reviewed labels for session %d: %s', ...
            review_session.session_index, review_session.session_label));
        drawnow;
        try
            review_session = synchronize_review_session_status(review_session);
            sync_current_session_cache();
            unverified_count = review_session.total_connections_identified - review_session.total_connections_verified;

            if unverified_count > 0
                choice = questdlg(sprintf('This session still has %d unverified candidate connections. Write labels anyway?', ...
                    unverified_count), ...
                    'Write Session Labels', 'Write Anyway', 'Cancel', 'Cancel');
                if ~strcmp(choice, 'Write Anyway')
                    return
                end
            end

            cell_classification = compute_reviewed_cell_labels(review_session);
            curated_session = build_curated_session(review_session, current_review_file);
            curated_file = fullfile(curated_session_root, sprintf('Session%03d-Animal%s-Day%s-Curated.mat', ...
                review_session.session_index, num2str(review_session.animal), num2str(review_session.day)));
            set(txt_status, 'String', sprintf(['Started writing reviewed labels for session %d: %s\n' ...
                'Exporting accepted connection plots...'], ...
                review_session.session_index, review_session.session_label));
            drawnow;
            writeback_export_options = struct('clearExisting', true);
            [accepted_plot_folder, accepted_plot_count] = export_accepted_connection_plots( ...
                review_session, accepted_plot_root, writeback_export_options);
            curated_session.accepted_plot_folder = accepted_plot_folder;
            curated_session.accepted_plot_count = accepted_plot_count;

            set(txt_status, 'String', sprintf(['Started writing reviewed labels for session %d: %s\n' ...
                'Saving curated session and updating All_Cells_combined classification...'], ...
                review_session.session_index, review_session.session_label));
            drawnow;
            save(curated_file, 'curated_session', '-v7.3');
            write_session_connection_exports(curated_session, curated_file);

            loaded_all = load(all_cells_file, 'All_Cells_combined');
            All_Cells_combined = loaded_all.All_Cells_combined;
            All_Cells_combined(review_session.session_index).CCGbased_CellClassfication = cell_classification;
            save(all_cells_file, 'All_Cells_combined', '-v7.3');

            review_session.reviewed_cell_classification = cell_classification;
            review_session.writeback_complete = true;
            review_session.last_writeback = char(datetime('now'));
            review_session.curated_file = curated_file;
            review_session.accepted_plot_folder = accepted_plot_folder;
            review_session.accepted_plot_count = accepted_plot_count;
            review_session = synchronize_review_session_status(review_session);
            sync_current_session_cache();
            session_dirty(current_file_idx) = true;
            save_current_session(true);

            set(txt_status, 'String', sprintf(['Started writing reviewed labels for session %d: %s\n' ...
                'Updating global curated summary...'], ...
                review_session.session_index, review_session.session_label));
            drawnow;
            update_global_curated_connectivity(curated_session_root, curated_root);
            global_summary = build_global_curated_summary(all_cells_file, curated_session_root);
            latest_summary_file = save_global_curated_summary(global_summary, summary_root);

            refresh_selectors();
            render_pair();
            set(txt_status, 'String', sprintf(['Wrote reviewed labels for session %d.\n' ...
                'Curated connectivity: %s\nAccepted plots: %s\nGlobal summary: %s'], ...
                review_session.session_index, curated_file, empty_as_none(accepted_plot_folder), latest_summary_file));
        catch exception
            set(txt_status, 'String', sprintf('Write Session Labels failed: %s', exception.message));
            errordlg(exception.message, 'Write Session Labels Failed');
        end
    end

    function write_global_summary_callback(~, ~)
        save_current_session(false);
        update_global_curated_connectivity(curated_session_root, curated_root);
        global_summary = build_global_curated_summary(all_cells_file, curated_session_root);
        latest_summary_file = save_global_curated_summary(global_summary, summary_root);
        refresh_selectors();
        render_pair();
        set(txt_status, 'String', sprintf('Updated global curated summary: %s', latest_summary_file));
    end

    function examine_presynaptic_cell(~, ~)
        if isempty(review_session.review_pairs)
            return
        end

        pair = review_session.review_pairs(current_idx);
        set(txt_status, 'String', sprintf('Preparing presynaptic cell view for %s...', ...
            review_session.tt_files{pair.pre}));
        drawnow;

        try
            presyn_fig = render_presynaptic_cell_view(review_session, pair.pre, presyn_fig);
            set(txt_status, 'String', sprintf('Opened presynaptic cell view for %s', ...
                review_session.tt_files{pair.pre}));
        catch exception
            set(txt_status, 'String', sprintf('Presynaptic cell view failed: %s', exception.message));
            errordlg(exception.message, 'Presynaptic Cell View');
        end
    end

    function close_reviewer(~, ~)
        save_current_session(false);
        delete(fig);
    end

    function load_session(target_idx, reset_to_first)
        if target_idx < 1 || target_idx > numel(review_files)
            return
        end

        save_current_session(false);

        [loaded_session, upgraded] = load_review_session_file(review_files{target_idx});
        review_session = loaded_session;
        current_review_file = review_files{target_idx};
        current_file_idx = target_idx;
        current_animal_idx = find_animal_index_for_file(session_cache, animal_groups, current_file_idx);
        active_file_indices = animal_groups{current_animal_idx};

        if upgraded
            save(current_review_file, 'review_session', '-v7.3');
        end

        sync_current_session_cache();

        if reset_to_first || isempty(review_session.review_pairs)
            current_idx = initial_pair_index(review_session, process_mode);
        else
            current_idx = min(max(current_idx, 1), max(numel(review_session.review_pairs), 1));
        end

        refresh_selectors();
        render_pair();
    end

    function save_current_session(force_save)
        if isempty(current_review_file)
            return
        end

        if force_save || session_dirty(current_file_idx)
            review_session = synchronize_review_session_status(review_session);
            save(current_review_file, 'review_session', '-v7.3');
            sync_current_session_cache();
            session_dirty(current_file_idx) = false;
        end
    end

    function refresh_selectors()
        sync_current_session_cache();
        session_strings = build_session_strings(session_cache, active_file_indices, current_file_idx);
        local_idx = current_active_session_position();
        set(popup_session, 'String', session_strings, 'Value', local_idx);

        pair_strings = build_pair_strings(review_session);
        set(popup_pair, 'String', pair_strings, 'Value', max(1, min(current_idx, numel(pair_strings))));

        has_pairs = ~isempty(review_session.review_pairs);
        pair_enable = on_off(has_pairs);
        set(btn_prev_pair, 'Enable', pair_enable);
        set(btn_next_pair, 'Enable', pair_enable);
        set(btn_accept, 'Enable', pair_enable);
        set(btn_discard, 'Enable', pair_enable);
        set(btn_switch, 'Enable', pair_enable);
        set(btn_set_exc, 'Enable', pair_enable);
        set(btn_set_inh, 'Enable', pair_enable);
        set(btn_examine_pre, 'Enable', pair_enable);
        set(popup_pair, 'Enable', pair_enable);

        set(btn_prev_session, 'Enable', on_off(local_idx > 1));
        set(btn_next_session, 'Enable', on_off(local_idx < numel(active_file_indices)));
        set(btn_next_animal, 'Enable', on_off(current_animal_idx < numel(animal_groups)));
        set(txt_mode, 'String', build_mode_string(process_mode, numel(review_files), ...
            skipped_mode_count, skipped_empty_count, animal_labels{current_animal_idx}, ...
            current_animal_idx, numel(animal_groups), numel(active_file_indices)));
    end

    function sync_current_session_cache()
        if isempty(current_review_file) || isempty(review_session)
            return
        end
        session_cache(current_file_idx) = make_review_session_cache_entry(current_review_file, review_session);
    end

    function local_idx = current_active_session_position()
        local_idx = find(active_file_indices == current_file_idx, 1, 'first');
        if isempty(local_idx)
            local_idx = 1;
        end
    end

    function switch_to_animal(target_animal_idx)
        if target_animal_idx < 1 || target_animal_idx > numel(animal_groups)
            return
        end
        current_animal_idx = target_animal_idx;
        active_file_indices = animal_groups{current_animal_idx};
        load_session(active_file_indices(1), true);
    end

    function render_pair()
        review_session = synchronize_review_session_status(review_session);
        sync_current_session_cache();
        cell_classification = compute_reviewed_cell_labels(review_session);

        if isempty(review_session.review_pairs)
            render_empty_axes(ax_full, 'Detector CCG (1.0 ms bins)');
            render_empty_axes(ax_zoom, 'Detector CCG (0.2 ms bins)');
            render_empty_axes(ax_raw, 'Raw CCG (0.1 ms, display-only)');

            info_lines = {
                sprintf('%s | No identified candidate connections in this session', review_session.session_label)
                sprintf('Pairs evaluated: %d | Identified: %d | Verified: %d | Writeback complete: %s', ...
                    review_session.total_pairs_evaluated, review_session.total_connections_identified, ...
                    review_session.total_connections_verified, ...
                    ternary(review_session.writeback_complete, 'yes', 'no'))
                sprintf('Cell labels now: excitatory %d | inhibitory %d | mixed %d | unlabeled %d', ...
                    sum(cell_classification == 1), sum(cell_classification == 2), ...
                    sum(cell_classification == 3), sum(cell_classification == 0))
                };

            set(txt_info, 'String', sprintf('%s\n%s\n%s', info_lines{1}, info_lines{2}, info_lines{3}));
            set(txt_status, 'String', sprintf('Viewing %s', current_review_file));
            return
        end

        pair = review_session.review_pairs(current_idx);
        plot_type = effective_plot_type(pair);
        plot_bundle = build_plot_bundle(review_session, pair, plot_type);

        render_ccg_axes(ax_full, plot_bundle.full_t_ms, plot_bundle.full_y, ...
            plot_bundle.full_pred_t_ms, plot_bundle.full_pred_y, ...
            plot_bundle.full_exc_upper_y, plot_bundle.full_inh_lower_y, ...
            plot_bundle.full_xlim_ms, plot_bundle.sig_window_exc_ms, plot_bundle.sig_window_inh_ms, ...
            sprintf('Detector CCG (%0.1f ms bins)', review_session.plot_config.plot_bin_full_s * 1000));

        render_ccg_axes(ax_zoom, plot_bundle.zoom_t_ms, plot_bundle.zoom_y, ...
            plot_bundle.zoom_pred_t_ms, plot_bundle.zoom_pred_y, ...
            plot_bundle.zoom_exc_upper_y, plot_bundle.zoom_inh_lower_y, ...
            plot_bundle.zoom_xlim_ms, plot_bundle.sig_window_exc_ms, plot_bundle.sig_window_inh_ms, ...
            sprintf('Detector CCG (%0.1f ms bins)', plot_bundle.zoom_bin_ms));

        render_ccg_axes(ax_raw, plot_bundle.raw_t_ms, plot_bundle.raw_y, ...
            [], [], [], [], ...
            plot_bundle.zoom_xlim_ms, plot_bundle.sig_window_exc_ms, plot_bundle.sig_window_inh_ms, ...
            'Raw CCG (0.1 ms, display-only)');

        local_session_idx = current_active_session_position();
        info_lines = {
            sprintf('%s | Animal session %d / %d | Pair %d / %d | Pre %s (FR %.3f Hz) -> Post %s (FR %.3f Hz)', ...
                review_session.session_label, local_session_idx, numel(active_file_indices), ...
                current_idx, numel(review_session.review_pairs), ...
                review_session.tt_files{pair.pre}, review_session.firing_rates(pair.pre), ...
                review_session.tt_files{pair.post}, review_session.firing_rates(pair.post))
            sprintf('Auto type: %s | Current decision: %s | Verified: %s | Pair scope: %s | Rebound override: %s | Last reviewed: %s', ...
                type_label(pair.auto_type), decision_label(pair.current_decision), ...
                ternary(pair.verified, 'yes', 'no'), ...
                ternary(pair.same_tt_pair, 'same-TT', 'cross-TT'), ...
                ternary(pair.rebound_override, 'yes', 'no'), ...
                empty_as_none(pair.last_reviewed))
            sprintf('Pairs evaluated: %d | Identified: %d | Verified: %d | Accepted excitatory: %d | Accepted inhibitory: %d | Cell labels now: excitatory %d | inhibitory %d | mixed %d | unlabeled %d', ...
                review_session.total_pairs_evaluated, review_session.total_connections_identified, ...
                review_session.total_connections_verified, ...
                sum([review_session.review_pairs.current_decision] == 1), ...
                sum([review_session.review_pairs.current_decision] == 2), ...
                sum(cell_classification == 1), sum(cell_classification == 2), ...
                sum(cell_classification == 3), sum(cell_classification == 0))
            };

        set(txt_info, 'String', sprintf('%s\n%s\n%s', info_lines{1}, info_lines{2}, info_lines{3}));
        set(txt_status, 'String', sprintf('Viewing %s', current_review_file));
    end
end

function path_options = normalize_path_options(export_options)

if isstruct(export_options) && numel(export_options) == 1
    path_options = export_options;
else
    path_options = struct();
end
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

function add_configured_dependency_paths(config, path_options)

cellExplorerPath = char(string(get_option_value(path_options, 'cellExplorerPath', config.cellExplorerPath)));
mclustPath = char(string(get_option_value(path_options, 'mclustPath', config.mclustPath)));

add_dependency_path(cellExplorerPath, 'CellExplorer');
add_dependency_path(mclustPath, 'MClust');

additionalPaths = get_option_value(path_options, 'additionalPaths', {});
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

function reviewRoot = default_review_root()

repoRoot = fullfile(fileparts(mfilename('fullpath')), '..', '..');
reviewRoot = fullfile(repoRoot, 'Data', 'CCG_review_sessions');
end

function curatedRoot = default_curated_root()

repoRoot = fullfile(fileparts(mfilename('fullpath')), '..', '..');
curatedRoot = fullfile(repoRoot, 'Data', 'CCG_curated_connections');
end

function cellPairsRoot = default_cell_pairs_root()

repoRoot = fullfile(fileparts(mfilename('fullpath')), '..', '..');
cellPairsRoot = fullfile(repoRoot, 'Results', 'CellPairs');
end

function allCellsPath = resolve_all_cells_path()

repoRoot = fullfile(fileparts(mfilename('fullpath')), '..', '..');
candidatePaths = { ...
    fullfile(repoRoot, 'All_Cells_combined.mat'), ...
    fullfile(repoRoot, 'Data', 'All_Cells_combined.mat'), ...
    fullfile(pwd, 'All_Cells_combined.mat')};

allCellsPath = first_existing_file(candidatePaths);
if isempty(allCellsPath)
    allCellsPath = fullfile(repoRoot, 'Data', 'All_Cells_combined.mat');
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
    sessionInfoPath = fullfile(repoRoot, 'Data', 'sessionInfo.mat');
end
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

function value = get_option_value(settings_struct, field_name, default_value)

value = default_value;
if nargin < 1 || isempty(settings_struct) || ~isstruct(settings_struct) || numel(settings_struct) ~= 1
    return
end

if isfield(settings_struct, field_name) && ~isempty(settings_struct.(field_name))
    value = settings_struct.(field_name);
end
end

function [review_files, process_mode] = select_review_targets(review_root)

review_files = {};
process_mode = '';

selection_choice = questdlg( ...
    'Load a whole review folder or choose individual review files?', ...
    'CCG Review Input', ...
    'Folder', 'Files', 'Cancel', 'Folder');

if isempty(selection_choice) || strcmp(selection_choice, 'Cancel')
    return
end

mode_choice = questdlg( ...
    'Process only new/unverified sessions, or re-open all selected sessions?', ...
    'CCG Review Mode', ...
    'Only New/Unverified', 'Reprocess All', 'Cancel', 'Only New/Unverified');

if isempty(mode_choice) || strcmp(mode_choice, 'Cancel')
    return
end

if strcmp(mode_choice, 'Reprocess All')
    process_mode = 'reprocess_all';
else
    process_mode = 'only_new';
end

if strcmp(selection_choice, 'Folder')
    selected_folder = uigetdir(review_root, 'Select a CCG review folder');
    if isequal(selected_folder, 0)
        review_files = {};
        return
    end
    review_files = resolve_review_files(selected_folder, review_root);
else
    [fname, fpath] = uigetfile(fullfile(review_root, '*.mat'), ...
        'Select one or more CCG review files', 'MultiSelect', 'on');
    if isequal(fname, 0)
        review_files = {};
        return
    end
    if iscell(fname)
        review_files = fullfile(fpath, fname);
    else
        review_files = {fullfile(fpath, fname)};
    end
end
end

function review_files = resolve_review_files(review_input, review_root)

if isempty(review_input)
    review_files = {};
    return
end

if isstring(review_input)
    review_input = cellstr(review_input);
end

if iscell(review_input)
    review_files = review_input(:);
elseif isfolder(review_input)
    listing = dir(fullfile(review_input, '*.mat'));
    review_files = fullfile({listing.folder}, {listing.name});
elseif exist(review_input, 'file')
    review_files = {review_input};
elseif exist(fullfile(review_root, review_input), 'file')
    review_files = {fullfile(review_root, review_input)};
else
    error('Could not resolve review input: %s', review_input);
end

review_files = sort(review_files(:));
end

function [filtered_files, session_cache, skipped_mode_count, skipped_empty_count, skipped_empty_summary] = ...
        filter_review_files(review_files, process_mode)

filtered_files = {};
session_cache = repmat(empty_review_session_cache_entry(), 0, 1);
skipped_mode_count = 0;
skipped_empty_count = 0;
skipped_empty_summary = struct('session_count', 0, 'pairs_evaluated', 0);

for idx = 1:numel(review_files)
    [review_session, upgraded] = load_review_session_file(review_files{idx});
    if upgraded
        save(review_files{idx}, 'review_session', '-v7.3');
    end

    if isempty(review_session.review_pairs)
        skipped_empty_count = skipped_empty_count + 1;
        skipped_empty_summary.session_count = skipped_empty_summary.session_count + 1;
        skipped_empty_summary.pairs_evaluated = skipped_empty_summary.pairs_evaluated + review_session.total_pairs_evaluated;
        continue
    end

    if session_requires_processing(review_session, process_mode)
        filtered_files{end+1,1} = review_files{idx}; %#ok<AGROW>
        session_cache(end+1,1) = make_review_session_cache_entry(review_files{idx}, review_session); %#ok<AGROW>
    else
        skipped_mode_count = skipped_mode_count + 1;
    end
end
end

function tf = session_requires_processing(review_session, process_mode)

review_session = synchronize_review_session_status(review_session);

switch process_mode
    case 'reprocess_all'
        tf = true;
    otherwise
        tf = ~review_session.review_complete || ~review_session.writeback_complete;
end
end

function batch_summary = export_accepted_connection_plots_from_review_files(review_files, accepted_plot_root, export_options)

export_options = normalize_accepted_plot_export_options(export_options);
ensure_dir(accepted_plot_root);

batch_summary = struct();
batch_summary.review_file_count = numel(review_files);
batch_summary.exported_session_count = 0;
batch_summary.skipped_empty_count = 0;
batch_summary.skipped_no_accepted_count = 0;
batch_summary.accepted_pair_count = 0;
batch_summary.output_root = accepted_plot_root;
batch_summary.include_original = export_options.includeOriginal;
batch_summary.clear_existing = export_options.clearExisting;
batch_summary.display_variants = export_options.displayVariants;
batch_summary.failures = repmat(struct('review_file', '', 'message', ''), 0, 1);

fprintf('\n=== Accepted CCG Plot Batch Export ===\n');
fprintf('Review files: %d\n', numel(review_files));
fprintf('Output root: %s\n', accepted_plot_root);
fprintf('Original detector plots: %s\n', ternary(export_options.includeOriginal, 'yes', 'no'));
fprintf('Clear existing session folders: %s\n', ternary(export_options.clearExisting, 'yes', 'no'));
for variant_idx = 1:numel(export_options.displayVariants)
    variant = export_options.displayVariants(variant_idx);
    fprintf('Display variant %d: %s | %.3g ms bins', ...
        variant_idx, variant.suffix, variant.zoomBinMs);
    if ~isempty(variant.zoomWindowMs)
        fprintf(' | window %.3g to %.3g ms', variant.zoomWindowMs(1), variant.zoomWindowMs(2));
    end
    fprintf('\n');
end

for file_idx = 1:numel(review_files)
    review_file = review_files{file_idx};
    try
        [review_session, ~] = load_review_session_file(review_file);
        review_session = synchronize_review_session_status(review_session);

        if isempty(review_session.review_pairs)
            batch_summary.skipped_empty_count = batch_summary.skipped_empty_count + 1;
            continue
        end

        accepted_pair_count = sum([review_session.review_pairs.current_decision] ~= 0);
        if accepted_pair_count == 0
            batch_summary.skipped_no_accepted_count = batch_summary.skipped_no_accepted_count + 1;
            continue
        end

        [accepted_plot_folder, accepted_plot_count] = export_accepted_connection_plots( ...
            review_session, accepted_plot_root, export_options);

        batch_summary.exported_session_count = batch_summary.exported_session_count + 1;
        batch_summary.accepted_pair_count = batch_summary.accepted_pair_count + accepted_plot_count;
        fprintf('  [%d/%d] Session %d: exported %d accepted pairs -> %s\n', ...
            file_idx, numel(review_files), review_session.session_index, ...
            accepted_plot_count, accepted_plot_folder);
    catch exception
        batch_summary.failures(end+1,1) = struct( ...
            'review_file', review_file, ...
            'message', exception.message);
        warning('Accepted plot export failed for %s: %s', review_file, exception.message);
    end
end

fprintf('Export complete: %d sessions, %d accepted pairs, %d failures.\n', ...
    batch_summary.exported_session_count, batch_summary.accepted_pair_count, ...
    numel(batch_summary.failures));
end

function idx = initial_pair_index(review_session, process_mode)

if isempty(review_session.review_pairs)
    idx = 1;
    return
end

if strcmp(process_mode, 'only_new')
    first_unverified = find(~[review_session.review_pairs.verified], 1, 'first');
    if ~isempty(first_unverified)
        idx = first_unverified;
        return
    end
end

idx = 1;
end

function [review_session, upgraded] = load_review_session_file(review_file)

loaded = load(review_file, 'review_session');
review_session = loaded.review_session;
[review_session, upgraded] = upgrade_review_session(review_session);
end

function [review_session, upgraded] = upgrade_review_session(review_session)

upgraded = false;

if ~isfield(review_session, 'version')
    review_session.version = 2;
    upgraded = true;
end

if ~isfield(review_session, 'review_pairs') || isempty(review_session.review_pairs)
    review_session.review_pairs = empty_review_pairs();
    upgraded = true;
end

for idx = 1:numel(review_session.review_pairs)
    if ~isfield(review_session.review_pairs, 'verified')
        review_session.review_pairs(idx).verified = false;
        upgraded = true;
    elseif isempty(review_session.review_pairs(idx).verified)
        review_session.review_pairs(idx).verified = false;
        upgraded = true;
    else
        review_session.review_pairs(idx).verified = logical(review_session.review_pairs(idx).verified);
    end

    if ~isfield(review_session.review_pairs, 'last_reviewed')
        review_session.review_pairs(idx).last_reviewed = '';
        upgraded = true;
    elseif isempty(review_session.review_pairs(idx).last_reviewed)
        review_session.review_pairs(idx).last_reviewed = '';
        upgraded = true;
    end

    if ~isfield(review_session.review_pairs, 'current_decision') || isempty(review_session.review_pairs(idx).current_decision)
        review_session.review_pairs(idx).current_decision = review_session.review_pairs(idx).initial_decision;
        upgraded = true;
    end
end

if ~isfield(review_session, 'candidate_signature')
    review_session.candidate_signature = make_candidate_signature(review_session.review_pairs);
    upgraded = true;
end
if ~isfield(review_session, 'settings_signature')
    review_session.settings_signature = '';
    upgraded = true;
end
if ~isfield(review_session, 'behavior_session_dirs')
    review_session.behavior_session_dirs = {};
    upgraded = true;
end
if ~isfield(review_session, 'behavior_spike_times')
    review_session.behavior_spike_times = {};
    upgraded = true;
end
if ~isfield(review_session, 'total_connections_identified')
    review_session.total_connections_identified = numel(review_session.review_pairs);
    upgraded = true;
end
if ~isfield(review_session, 'total_pairs_evaluated')
    review_session.total_pairs_evaluated = resolve_total_pairs_evaluated(review_session);
    upgraded = true;
end
if ~isfield(review_session, 'total_connections_tested')
    review_session.total_connections_tested = numel(review_session.review_pairs);
    upgraded = true;
end
if ~isfield(review_session, 'total_connections_verified')
    review_session.total_connections_verified = sum([review_session.review_pairs.verified]);
    upgraded = true;
end
if ~isfield(review_session, 'review_complete')
    review_session.review_complete = review_session.total_connections_verified == review_session.total_connections_identified;
    upgraded = true;
end
if ~isfield(review_session, 'writeback_complete')
    review_session.writeback_complete = isfield(review_session, 'last_writeback') && ~isempty(review_session.last_writeback);
    upgraded = true;
end
if ~isfield(review_session, 'last_writeback')
    review_session.last_writeback = '';
    upgraded = true;
end
if ~isfield(review_session, 'curated_file')
    review_session.curated_file = '';
    upgraded = true;
end
if ~isfield(review_session, 'reviewed_cell_classification')
    review_session.reviewed_cell_classification = [];
    upgraded = true;
end

review_session = synchronize_review_session_status(review_session);
end

function review_session = synchronize_review_session_status(review_session)

if ~isfield(review_session, 'review_pairs') || isempty(review_session.review_pairs)
    review_session.review_pairs = empty_review_pairs();
end

review_session.total_connections_identified = numel(review_session.review_pairs);
review_session.total_connections_tested = review_session.total_connections_identified;
review_session.total_pairs_evaluated = resolve_total_pairs_evaluated(review_session);
if isempty(review_session.review_pairs)
    review_session.total_connections_verified = 0;
else
    review_session.total_connections_verified = sum([review_session.review_pairs.verified]);
end
review_session.review_complete = review_session.total_connections_verified == review_session.total_connections_identified;
review_session.candidate_signature = make_candidate_signature(review_session.review_pairs);

if ~isfield(review_session, 'writeback_complete') || isempty(review_session.writeback_complete)
    review_session.writeback_complete = false;
end
if ~isfield(review_session, 'last_writeback')
    review_session.last_writeback = '';
end
if ~isfield(review_session, 'curated_file')
    review_session.curated_file = '';
end
if ~isfield(review_session, 'settings_signature')
    review_session.settings_signature = '';
end
if ~isfield(review_session, 'behavior_session_dirs')
    review_session.behavior_session_dirs = {};
end
if ~isfield(review_session, 'behavior_spike_times')
    review_session.behavior_spike_times = {};
end
end

function review_pairs = empty_review_pairs()

review_pairs = repmat(struct( ...
    'pre', 0, ...
    'post', 0, ...
    'auto_type', 0, ...
    'initial_decision', 0, ...
    'current_decision', 0, ...
    'same_tt_pair', false, ...
    'rebound_override', false, ...
    'verified', false, ...
    'last_reviewed', '' ...
    ), 0, 1);
end

function strings = build_session_strings(session_cache, active_file_indices, current_file_idx)

strings = cell(numel(active_file_indices), 1);

for idx = 1:numel(active_file_indices)
    cache_entry = session_cache(active_file_indices(idx));

    strings{idx} = sprintf('%s | identified %d | verified %d | evaluated %d | writeback %s', ...
        cache_entry.session_label, cache_entry.identified_connections, cache_entry.verified_connections, ...
        cache_entry.total_pairs_evaluated, ternary(cache_entry.writeback_complete, 'yes', 'no'));

    if active_file_indices(idx) == current_file_idx
        strings{idx} = [strings{idx} ' | current'];
    end
end
end

function cache_entry = empty_review_session_cache_entry()

cache_entry = struct( ...
    'review_file', '', ...
    'session_index', 0, ...
    'session_label', '', ...
    'animal', '', ...
    'animal_key', '', ...
    'animal_label', '', ...
    'identified_connections', 0, ...
    'verified_connections', 0, ...
    'total_pairs_evaluated', 0, ...
    'writeback_complete', false, ...
    'review_complete', false ...
    );
end

function cache_entry = make_review_session_cache_entry(review_file, review_session)

review_session = synchronize_review_session_status(review_session);
animal_label = stringify_animal_value(review_session.animal);

cache_entry = empty_review_session_cache_entry();
cache_entry.review_file = review_file;
cache_entry.session_index = review_session.session_index;
cache_entry.session_label = review_session.session_label;
cache_entry.animal = review_session.animal;
cache_entry.animal_key = animal_label;
cache_entry.animal_label = animal_label;
cache_entry.identified_connections = review_session.total_connections_identified;
cache_entry.verified_connections = review_session.total_connections_verified;
cache_entry.total_pairs_evaluated = review_session.total_pairs_evaluated;
cache_entry.writeback_complete = review_session.writeback_complete;
cache_entry.review_complete = review_session.review_complete;
end

function animal_groups = build_animal_groups(session_cache)

animal_groups = {};
if isempty(session_cache)
    return
end

animal_keys = {session_cache.animal_key};
[unique_keys, ~, key_idx] = unique(animal_keys, 'stable');

animal_groups = cell(numel(unique_keys), 1);
for idx = 1:numel(unique_keys)
    animal_groups{idx} = find(key_idx == idx);
end
end

function animal_labels = build_animal_labels(session_cache, animal_groups)

animal_labels = cell(numel(animal_groups), 1);
for idx = 1:numel(animal_groups)
    first_idx = animal_groups{idx}(1);
    animal_labels{idx} = session_cache(first_idx).animal_label;
end
end

function animal_idx = find_animal_index_for_file(session_cache, animal_groups, file_idx)

animal_idx = 1;
for idx = 1:numel(animal_groups)
    if any(animal_groups{idx} == file_idx)
        animal_idx = idx;
        return
    end
end

if ~isempty(session_cache)
    file_key = session_cache(file_idx).animal_key;
    all_labels = build_animal_labels(session_cache, animal_groups);
    animal_idx = find(strcmp(all_labels, file_key), 1, 'first');
    if isempty(animal_idx)
        animal_idx = 1;
    end
end
end

function label_text = build_mode_string(process_mode, file_count, skipped_mode_count, skipped_empty_count, ...
        animal_label, current_animal_idx, total_animals, animal_session_count)

label_text = sprintf(['Mode: %s | Sessions loaded: %d | Skipped by mode: %d | ' ...
    'Skipped empty: %d | Animal %d/%d: %s (%d sessions)'], ...
    process_mode_label(process_mode), file_count, skipped_mode_count, skipped_empty_count, ...
    current_animal_idx, total_animals, animal_label, animal_session_count);
end

function value = resolve_total_pairs_evaluated(review_session)

value = [];

if isfield(review_session, 'auto_pair_counts') && ~isempty(review_session.auto_pair_counts)
    counts = review_session.auto_pair_counts;
    if isfield(counts, 'total_pairs_evaluated') && isfinite(counts.total_pairs_evaluated)
        value = double(counts.total_pairs_evaluated);
    elseif isfield(counts, 'total_connections_tested') && ...
            isfield(counts, 'identified_candidate_pairs') && ...
            isfinite(counts.total_connections_tested)
        value = double(counts.total_connections_tested);
    end
end

if isempty(value)
    tested_unit_count = resolve_tested_unit_count(review_session);
    value = tested_unit_count * max(tested_unit_count - 1, 0);
end
end

function tested_unit_count = resolve_tested_unit_count(review_session)

tested_unit_count = 0;

if isfield(review_session, 'spike_times') && iscell(review_session.spike_times) && ~isempty(review_session.spike_times)
    tested_unit_count = sum(cellfun(@(ts) ~isempty(ts) && any(isfinite(ts(:))), review_session.spike_times));
end

if tested_unit_count == 0 && isfield(review_session, 'n_units') && ~isempty(review_session.n_units)
    tested_unit_count = double(review_session.n_units);
end
end

function text_value = stringify_animal_value(raw_value)

if isnumeric(raw_value)
    if isscalar(raw_value)
        text_value = num2str(raw_value);
    else
        text_value = mat2str(raw_value);
    end
elseif isstring(raw_value)
    text_value = char(raw_value);
else
    text_value = char(string(raw_value));
end
end

function strings = build_pair_strings(review_session)

if isempty(review_session.review_pairs)
    strings = {'No identified candidate connections'};
    return
end

strings = cell(numel(review_session.review_pairs), 1);
for idx = 1:numel(review_session.review_pairs)
    pair = review_session.review_pairs(idx);
    strings{idx} = sprintf('[%s] %02d | %s -> %s | auto %s | current %s | %s', ...
        ternary(pair.verified, 'V', ' '), idx, ...
        review_session.tt_files{pair.pre}, review_session.tt_files{pair.post}, ...
        type_label(pair.auto_type), decision_label(pair.current_decision), ...
        ternary(pair.same_tt_pair, 'same-TT', 'cross-TT'));
end
end

function plot_bundle = build_plot_bundle(review_session, pair, plot_type, display_zoom_bin_ms, display_zoom_window_ms)

if nargin < 4 || isempty(display_zoom_bin_ms)
    display_zoom_bin_ms = min(review_session.plot_config.plot_bin_zoom_exc_s, ...
        review_session.plot_config.plot_bin_zoom_inh_s) * 1000;
end
if nargin < 5 || isempty(display_zoom_window_ms)
    display_zoom_window_ms = resolve_common_zoom_window(review_session.plot_config);
end

if plot_type == 2 && ~isempty(review_session.mono_res_inhibitory)
    primary_plot = review_session.mono_res_inhibitory;
elseif ~isempty(review_session.mono_res_excitatory)
    primary_plot = review_session.mono_res_excitatory;
else
    primary_plot = review_session.mono_res_inhibitory;
end

exc_plot = review_session.mono_res_excitatory;
inh_plot = review_session.mono_res_inhibitory;
zoom_window_ms = display_zoom_window_ms;

ccg_ts_ms = primary_plot.ccgTs(:) * 1000;
pred_ts_ms = primary_plot.ccgTs(:) * 1000;

[full_t_ms, full_y] = rebin_ccg_counts(ccg_ts_ms, primary_plot.ccgR(:,pair.pre,pair.post), ...
    review_session.plot_config.plot_bin_full_s * 1000, review_session.plot_config.full_window_ms);
[full_pred_t_ms, full_pred_y] = rebin_ccg_counts(pred_ts_ms, primary_plot.Pred(:,pair.pre,pair.post), ...
    review_session.plot_config.plot_bin_full_s * 1000, review_session.plot_config.full_window_ms);
[~, full_exc_upper_y] = rebin_ccg_counts(pred_ts_ms, exc_plot.Bounds(:,pair.pre,pair.post,1), ...
    review_session.plot_config.plot_bin_full_s * 1000, review_session.plot_config.full_window_ms);
[~, full_inh_lower_y] = rebin_ccg_counts(pred_ts_ms, inh_plot.Bounds(:,pair.pre,pair.post,2), ...
    review_session.plot_config.plot_bin_full_s * 1000, review_session.plot_config.full_window_ms);

[zoom_t_ms, zoom_y] = rebin_ccg_counts(ccg_ts_ms, primary_plot.ccgR(:,pair.pre,pair.post), ...
    display_zoom_bin_ms, zoom_window_ms);
[zoom_pred_t_ms, zoom_pred_y] = rebin_ccg_counts(pred_ts_ms, primary_plot.Pred(:,pair.pre,pair.post), ...
    display_zoom_bin_ms, zoom_window_ms);
[~, zoom_exc_upper_y] = rebin_ccg_counts(pred_ts_ms, exc_plot.Bounds(:,pair.pre,pair.post,1), ...
    display_zoom_bin_ms, zoom_window_ms);
[~, zoom_inh_lower_y] = rebin_ccg_counts(pred_ts_ms, inh_plot.Bounds(:,pair.pre,pair.post,2), ...
    display_zoom_bin_ms, zoom_window_ms);

[raw_y, raw_t_s] = CrossCorrel(review_session.spike_times{pair.pre}, review_session.spike_times{pair.post}, ...
    review_session.plot_config.plot_bin_raw_s, zoom_window_ms / 1000);

plot_bundle = struct();
plot_bundle.full_t_ms = full_t_ms;
plot_bundle.full_y = full_y;
plot_bundle.full_pred_t_ms = full_pred_t_ms;
plot_bundle.full_pred_y = full_pred_y;
plot_bundle.full_exc_upper_y = full_exc_upper_y;
plot_bundle.full_inh_lower_y = full_inh_lower_y;
plot_bundle.full_xlim_ms = review_session.plot_config.full_window_ms;
plot_bundle.zoom_t_ms = zoom_t_ms;
plot_bundle.zoom_y = zoom_y;
plot_bundle.zoom_pred_t_ms = zoom_pred_t_ms;
plot_bundle.zoom_pred_y = zoom_pred_y;
plot_bundle.zoom_exc_upper_y = zoom_exc_upper_y;
plot_bundle.zoom_inh_lower_y = zoom_inh_lower_y;
plot_bundle.zoom_xlim_ms = zoom_window_ms;
plot_bundle.zoom_bin_ms = display_zoom_bin_ms;
plot_bundle.raw_t_ms = raw_t_s(:) * 1000;
plot_bundle.raw_y = double(raw_y(:));
plot_bundle.sig_window_exc_ms = review_session.plot_config.sig_window_exc_ms;
plot_bundle.sig_window_inh_ms = review_session.plot_config.sig_window_inh_ms;
end

function render_ccg_axes(ax, bar_t_ms, bar_y, pred_t_ms, pred_y, exc_upper_y, inh_lower_y, ...
        xwin_ms, sig_window_exc_ms, sig_window_inh_ms, ax_title)

cla(ax)
bar(ax, bar_t_ms, bar_y, 1, 'FaceColor', [0.3 0.3 0.3], 'EdgeColor', 'none')
hold(ax, 'on')

yl = ylim(ax);
patch(ax, [sig_window_exc_ms(1) sig_window_exc_ms(2) sig_window_exc_ms(2) sig_window_exc_ms(1)], ...
    [yl(1) yl(1) yl(2) yl(2)], [0.2 0.4 0.9], 'FaceAlpha', 0.08, 'EdgeColor', 'none');
patch(ax, [sig_window_inh_ms(1) sig_window_inh_ms(2) sig_window_inh_ms(2) sig_window_inh_ms(1)], ...
    [yl(1) yl(1) yl(2) yl(2)], [0.9 0.3 0.3], 'FaceAlpha', 0.06, 'EdgeColor', 'none');

if ~isempty(pred_t_ms)
    plot(ax, pred_t_ms, pred_y, 'k-', 'LineWidth', 1.2)
    if ~isempty(exc_upper_y)
        plot(ax, pred_t_ms, exc_upper_y, '--', 'Color', [0.15 0.35 0.85], 'LineWidth', 0.9)
    end
    if ~isempty(inh_lower_y)
        plot(ax, pred_t_ms, inh_lower_y, '--', 'Color', [0.85 0.2 0.2], 'LineWidth', 0.9)
    end
end

xlim(ax, xwin_ms)
xlabel(ax, 'time [ms]');
ylabel(ax, 'Count');
title(ax, ax_title);
hold(ax, 'off')
end

function [accepted_plot_folder, accepted_plot_count] = export_accepted_connection_plots(review_session, accepted_plot_root, export_options)

if nargin < 3
    export_options = [];
end

review_session = synchronize_review_session_status(review_session);
export_options = normalize_accepted_plot_export_options(export_options);
accepted_plot_folder = '';
accepted_plot_count = 0;

if isempty(review_session.review_pairs)
    return
end

accepted_idx = find([review_session.review_pairs.current_decision] ~= 0);
session_plot_folder = fullfile(accepted_plot_root, build_accepted_plot_session_folder_name(review_session));

if export_options.clearExisting && exist(session_plot_folder, 'dir')
    rmdir(session_plot_folder, 's');
end

if isempty(accepted_idx)
    return
end

ensure_dir(session_plot_folder);
crossTT_folder = fullfile(session_plot_folder, 'CrossTT');
sameTT_folder = fullfile(session_plot_folder, 'SameTT');
ensure_dir(crossTT_folder);
ensure_dir(sameTT_folder);

for idx = 1:numel(accepted_idx)
    pair_idx = accepted_idx(idx);
    pair = review_session.review_pairs(pair_idx);
    plot_type = effective_plot_type(pair);
    plot_bundle = build_plot_bundle(review_session, pair, plot_type);

    outfolder_this = ternary(pair.same_tt_pair, sameTT_folder, crossTT_folder);
    current_type = sanitize_filename_component(type_label(pair.current_decision));
    pair_scope = ternary(pair.same_tt_pair, 'sameTT', 'crossTT');
    base_name = sprintf('Pair%03d-%s-to-%s-%s-%s', ...
        pair_idx, ...
        sanitize_filename_component(review_session.tt_files{pair.pre}), ...
        sanitize_filename_component(review_session.tt_files{pair.post}), ...
        current_type, pair_scope);

    if export_options.includeOriginal
        pdf_file = fullfile(outfolder_this, [base_name, '.pdf']);
        png_file = fullfile(outfolder_this, [base_name, '.png']);
        write_accepted_pair_plot_figure(review_session, pair, pair_idx, plot_bundle, ...
            pdf_file, png_file, 'Accepted CCG Pair Export');
    end

    for variant_idx = 1:numel(export_options.displayVariants)
        variant = export_options.displayVariants(variant_idx);
        plot_bundle_extra = build_plot_bundle(review_session, pair, plot_type, ...
            variant.zoomBinMs, variant.zoomWindowMs);
        extra_base_name = [base_name, '-', sanitize_filename_component(variant.suffix)];
        extra_pdf_file = fullfile(outfolder_this, [extra_base_name, '.pdf']);
        extra_png_file = fullfile(outfolder_this, [extra_base_name, '.png']);
        write_accepted_pair_plot_figure(review_session, pair, pair_idx, plot_bundle_extra, ...
            extra_pdf_file, extra_png_file, sprintf('Accepted CCG Pair Export %0.1f ms', variant.zoomBinMs));
    end

    accepted_plot_count = accepted_plot_count + 1;
end

accepted_plot_folder = session_plot_folder;
end

function export_options = normalize_accepted_plot_export_options(export_options)

if nargin < 1 || isempty(export_options)
    export_options = struct();
end

if isstruct(export_options) && isfield(export_options, 'zoomBinMs') && ...
        ~isfield(export_options, 'displayVariants')
    export_options = struct('includeOriginal', true, 'displayVariants', export_options);
elseif isstruct(export_options) && numel(export_options) > 1 && ...
        ~isfield(export_options, 'displayVariants')
    export_options = struct('includeOriginal', true, 'displayVariants', export_options(:));
end

if ~isfield(export_options, 'includeOriginal') || isempty(export_options.includeOriginal)
    export_options.includeOriginal = true;
else
    export_options.includeOriginal = logical(export_options.includeOriginal);
end

if ~isfield(export_options, 'clearExisting') || isempty(export_options.clearExisting)
    export_options.clearExisting = false;
else
    export_options.clearExisting = logical(export_options.clearExisting);
end

if ~isfield(export_options, 'displayVariants') || isempty(export_options.displayVariants)
    export_options.displayVariants = default_accepted_plot_display_variants();
else
    export_options.displayVariants = standardize_accepted_plot_display_variants(export_options.displayVariants);
end
end

function variants = default_accepted_plot_display_variants()

variants = struct( ...
    'suffix', {'Display0400us'}, ...
    'zoomBinMs', {0.4}, ...
    'zoomWindowMs', {[]} ...
    );
end

function variants = standardize_accepted_plot_display_variants(variants)

if ~isstruct(variants)
    error('Accepted plot display variants must be provided as a struct array.');
end

variants = variants(:);
for idx = 1:numel(variants)
    if ~isfield(variants, 'zoomBinMs')
        error('Each accepted plot display variant needs a finite zoomBinMs value.');
    end

    zoom_bin_ms = variants(idx).zoomBinMs;
    if isempty(zoom_bin_ms) || numel(zoom_bin_ms) ~= 1 || ~isfinite(zoom_bin_ms)
        error('Each accepted plot display variant needs a finite zoomBinMs value.');
    end
    variants(idx).zoomBinMs = double(zoom_bin_ms);

    if ~isfield(variants, 'zoomWindowMs')
        variants(idx).zoomWindowMs = [];
    end
    if ~isempty(variants(idx).zoomWindowMs) && ...
            (numel(variants(idx).zoomWindowMs) ~= 2 || any(~isfinite(variants(idx).zoomWindowMs)))
        error('zoomWindowMs must be empty or a two-element [min max] vector.');
    end
    variants(idx).zoomWindowMs = double(variants(idx).zoomWindowMs(:)');
    if ~isempty(variants(idx).zoomWindowMs) && variants(idx).zoomWindowMs(1) >= variants(idx).zoomWindowMs(2)
        error('zoomWindowMs must be ordered as [min max].');
    end

    if ~isfield(variants, 'suffix') || isempty(variants(idx).suffix)
        variants(idx).suffix = build_display_variant_suffix( ...
            variants(idx).zoomBinMs, variants(idx).zoomWindowMs);
    end
end
end

function suffix = build_display_variant_suffix(zoom_bin_ms, zoom_window_ms)

suffix = sprintf('Display%04dus', round(zoom_bin_ms * 1000));
if ~isempty(zoom_window_ms)
    suffix = sprintf('%s-Win%sTo%sms', suffix, ...
        signed_number_token(zoom_window_ms(1)), signed_number_token(zoom_window_ms(2)));
end
end

function token = signed_number_token(value)

token = sprintf('%+g', value);
token = strrep(token, '+', 'p');
token = strrep(token, '-', 'm');
token = strrep(token, '.', 'p');
end

function write_accepted_pair_plot_figure(review_session, pair, pair_idx, plot_bundle, pdf_file, png_file, fig_name)

fig = figure( ...
    'Visible', 'off', ...
    'Color', 'w', ...
    'Units', 'pixels', ...
    'Position', [120 120 1700 760], ...
    'PaperPositionMode', 'auto', ...
    'NumberTitle', 'off', ...
    'Name', fig_name ...
    );

try
    annotation(fig, 'textbox', [0.03 0.90 0.94 0.08], ...
        'String', sprintf(['%s | Pair %d | Pre %s (FR %.3f Hz) -> Post %s (FR %.3f Hz)\n' ...
        'Current %s | Auto %s | Pair scope %s | Rebound override %s\n' ...
        'Blue dashed = excitatory upper bound | Red dashed = inhibitory lower bound'], ...
        review_session.session_label, pair_idx, ...
        review_session.tt_files{pair.pre}, review_session.firing_rates(pair.pre), ...
        review_session.tt_files{pair.post}, review_session.firing_rates(pair.post), ...
        decision_label(pair.current_decision), type_label(pair.auto_type), ...
        ternary(pair.same_tt_pair, 'same-TT', 'cross-TT'), ...
        ternary(pair.rebound_override, 'yes', 'no')), ...
        'Interpreter', 'none', ...
        'EdgeColor', 'none', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'middle', ...
        'FontSize', 11, ...
        'Color', [0.15 0.15 0.15]);

    ax_full = axes('Parent', fig, 'Units', 'normalized', 'Position', [0.05 0.12 0.27 0.70]);
    ax_zoom = axes('Parent', fig, 'Units', 'normalized', 'Position', [0.37 0.12 0.27 0.70]);
    ax_raw = axes('Parent', fig, 'Units', 'normalized', 'Position', [0.69 0.12 0.27 0.70]);

    render_ccg_axes(ax_full, plot_bundle.full_t_ms, plot_bundle.full_y, ...
        plot_bundle.full_pred_t_ms, plot_bundle.full_pred_y, ...
        plot_bundle.full_exc_upper_y, plot_bundle.full_inh_lower_y, ...
        plot_bundle.full_xlim_ms, plot_bundle.sig_window_exc_ms, plot_bundle.sig_window_inh_ms, ...
        sprintf('Detector CCG (%0.1f ms bins)', review_session.plot_config.plot_bin_full_s * 1000));

    render_ccg_axes(ax_zoom, plot_bundle.zoom_t_ms, plot_bundle.zoom_y, ...
        plot_bundle.zoom_pred_t_ms, plot_bundle.zoom_pred_y, ...
        plot_bundle.zoom_exc_upper_y, plot_bundle.zoom_inh_lower_y, ...
        plot_bundle.zoom_xlim_ms, plot_bundle.sig_window_exc_ms, plot_bundle.sig_window_inh_ms, ...
        sprintf('Detector CCG (%0.1f ms bins)', plot_bundle.zoom_bin_ms));

    render_ccg_axes(ax_raw, plot_bundle.raw_t_ms, plot_bundle.raw_y, ...
        [], [], [], [], ...
        plot_bundle.zoom_xlim_ms, plot_bundle.sig_window_exc_ms, plot_bundle.sig_window_inh_ms, ...
        'Raw CCG (0.1 ms, display-only)');

    exportgraphics(fig, pdf_file, 'ContentType', 'vector');
    exportgraphics(fig, png_file, 'Resolution', 300);
catch exception
    if ishandle(fig)
        close(fig);
    end
    rethrow(exception)
end

if ishandle(fig)
    close(fig);
end
end

function folder_name = build_accepted_plot_session_folder_name(review_session)

animal_text = sanitize_filename_component(stringify_animal_value(review_session.animal));
day_text = sanitize_filename_component(stringify_animal_value(review_session.day));
folder_name = sprintf('Session%03d-Animal%s-Day%s', review_session.session_index, animal_text, day_text);
end

function safe_text = sanitize_filename_component(raw_text)

safe_text = char(string(raw_text));
safe_text = regexprep(safe_text, '\s+', '_');
safe_text = regexprep(safe_text, '[^\w-]+', '_');
safe_text = regexprep(safe_text, '_+', '_');
safe_text = regexprep(safe_text, '^_+|_+$', '');

if isempty(safe_text)
    safe_text = 'item';
end
end

function zoom_window_ms = resolve_common_zoom_window(plot_config)

zoom_window_ms = [ ...
    min(plot_config.zoom_window_exc_ms(1), plot_config.zoom_window_inh_ms(1)) ...
    max(plot_config.zoom_window_exc_ms(2), plot_config.zoom_window_inh_ms(2)) ...
    ];
end

function render_empty_axes(ax, ax_title)

cla(ax)
title(ax, ax_title);
xlabel(ax, 'time [ms]');
ylabel(ax, 'Count');
text(ax, 0.5, 0.5, 'No identified candidate pairs', 'Units', 'normalized', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end

function cell_classification = compute_reviewed_cell_labels(review_session)

cell_classification = zeros(review_session.n_units, 1);

if isempty(review_session.review_pairs)
    return
end

accepted_exc = review_session.review_pairs([review_session.review_pairs.current_decision] == 1);
accepted_inh = review_session.review_pairs([review_session.review_pairs.current_decision] == 2);

exc_cells = [];
inh_cells = [];

if ~isempty(accepted_exc)
    exc_cells = [accepted_exc.pre]';
end
if ~isempty(accepted_inh)
    inh_cells = [accepted_inh.pre]';
end

all_cells = unique([exc_cells; inh_cells]);

for idx = 1:numel(all_cells)
    cell_id = all_cells(idx);
    exc_count = sum(exc_cells == cell_id);
    inh_count = sum(inh_cells == cell_id);

    if exc_count > 0 && inh_count == 0
        cell_classification(cell_id) = 1;
    elseif inh_count > 0 && exc_count == 0
        cell_classification(cell_id) = 2;
    else
        cell_classification(cell_id) = 3;
    end
end
end

function curated_session = build_curated_session(review_session, review_file)

review_session = synchronize_review_session_status(review_session);
cell_classification = compute_reviewed_cell_labels(review_session);
connectivity = build_reviewed_connectivity(review_session, cell_classification);

curated_session = struct();
curated_session.version = 1;
curated_session.session_index = review_session.session_index;
curated_session.session_label = review_session.session_label;
curated_session.mainDir = review_session.mainDir;
curated_session.animal = review_session.animal;
curated_session.day = review_session.day;
curated_session.n_units = review_session.n_units;
curated_session.tt_files = review_session.tt_files;
curated_session.firing_rates = review_session.firing_rates;
curated_session.review_file = review_file;
curated_session.curated_on = char(datetime('now'));
curated_session.settings_signature = review_session.settings_signature;
curated_session.total_connections_tested = review_session.total_pairs_evaluated;
curated_session.total_pairs_evaluated = review_session.total_pairs_evaluated;
curated_session.total_connections_identified = review_session.total_connections_identified;
curated_session.total_connections_verified = review_session.total_connections_verified;
curated_session.review_complete = review_session.review_complete;
curated_session.final_cell_classification = cell_classification;
curated_session.connectivity = connectivity;
curated_session.accepted_connections = connectivity.accepted_connections;
curated_session.accepted_connection_count = numel(connectivity.accepted_connections);
curated_session.accepted_plot_folder = '';
curated_session.accepted_plot_count = 0;
end

function connectivity = build_reviewed_connectivity(review_session, cell_classification)

accepted_pairs = review_session.review_pairs([review_session.review_pairs.current_decision] ~= 0);

accepted_connections = repmat(accepted_connection_template(), 0, 1);
for idx = 1:numel(accepted_pairs)
    pair = accepted_pairs(idx);
    accepted_connections(end+1,1) = struct( ... %#ok<AGROW>
        'session_index', review_session.session_index, ...
        'session_label', review_session.session_label, ...
        'animal', review_session.animal, ...
        'day', review_session.day, ...
        'pre', pair.pre, ...
        'post', pair.post, ...
        'pre_tt', review_session.tt_files{pair.pre}, ...
        'post_tt', review_session.tt_files{pair.post}, ...
        'pre_firing_rate', review_session.firing_rates(pair.pre), ...
        'post_firing_rate', review_session.firing_rates(pair.post), ...
        'connection_type_id', pair.current_decision, ...
        'connection_type', type_label(pair.current_decision), ...
        'auto_type_id', pair.auto_type, ...
        'auto_type', type_label(pair.auto_type), ...
        'initial_decision_id', pair.initial_decision, ...
        'initial_decision', decision_label(pair.initial_decision), ...
        'switched_from_auto', pair.current_decision ~= pair.initial_decision, ...
        'same_tt_pair', pair.same_tt_pair, ...
        'rebound_override', pair.rebound_override, ...
        'verified', pair.verified, ...
        'last_reviewed', pair.last_reviewed, ...
        'presynaptic_cell_label', cell_classification(pair.pre) ...
        );
end

accepted_exc_pairs_all = zeros(0,2);
accepted_inh_pairs_all = zeros(0,2);

if ~isempty(accepted_pairs)
    accepted_pre = [accepted_pairs.pre]';
    accepted_post = [accepted_pairs.post]';
    accepted_decision = [accepted_pairs.current_decision]';
    accepted_pre_label = cell_classification(accepted_pre);

    exc_mask = accepted_decision == 1 & accepted_pre_label == 1;
    inh_mask = accepted_decision == 2 & accepted_pre_label == 2;

    if any(exc_mask)
        accepted_exc_pairs_all = [accepted_pre(exc_mask) accepted_post(exc_mask)];
    end
    if any(inh_mask)
        accepted_inh_pairs_all = [accepted_pre(inh_mask) accepted_post(inh_mask)];
    end
end

[accepted_exc_pairs_crossTT, accepted_exc_pairs_sameTT] = split_pairs_by_tetrode(accepted_exc_pairs_all, review_session);
[accepted_inh_pairs_crossTT, accepted_inh_pairs_sameTT] = split_pairs_by_tetrode(accepted_inh_pairs_all, review_session);

conn_matrix = zeros(review_session.n_units);
for idx = 1:size(accepted_exc_pairs_crossTT, 1)
    conn_matrix(accepted_exc_pairs_crossTT(idx,2), accepted_exc_pairs_crossTT(idx,1)) = 1;
end
for idx = 1:size(accepted_inh_pairs_crossTT, 1)
    conn_matrix(accepted_inh_pairs_crossTT(idx,2), accepted_inh_pairs_crossTT(idx,1)) = -1;
end

connectivity = struct();
connectivity.accepted_connections = accepted_connections;
connectivity.accepted_exc_pairs_all = accepted_exc_pairs_all;
connectivity.accepted_inh_pairs_all = accepted_inh_pairs_all;
connectivity.accepted_exc_pairs_crossTT = accepted_exc_pairs_crossTT;
connectivity.accepted_exc_pairs_sameTT = accepted_exc_pairs_sameTT;
connectivity.accepted_inh_pairs_crossTT = accepted_inh_pairs_crossTT;
connectivity.accepted_inh_pairs_sameTT = accepted_inh_pairs_sameTT;
connectivity.crossTT_connectivity_matrix = conn_matrix;
connectivity.crossTT_connectivity_matrix_plot = conn_matrix + 2;
end

function [crossTT_pairs, sameTT_pairs] = split_pairs_by_tetrode(pairs, review_session)

if isempty(pairs)
    crossTT_pairs = zeros(0,2);
    sameTT_pairs = zeros(0,2);
    return
end

same_mask = false(size(pairs,1), 1);
for idx = 1:size(pairs,1)
    match_idx = find(arrayfun(@(p) p.pre == pairs(idx,1) && p.post == pairs(idx,2), review_session.review_pairs), 1, 'first');
    if ~isempty(match_idx)
        same_mask(idx) = review_session.review_pairs(match_idx).same_tt_pair;
    end
end

sameTT_pairs = pairs(same_mask, :);
crossTT_pairs = pairs(~same_mask, :);
end

function template = accepted_connection_template()

template = struct( ...
    'session_index', 0, ...
    'session_label', '', ...
    'animal', 0, ...
    'day', 0, ...
    'pre', 0, ...
    'post', 0, ...
    'pre_tt', '', ...
    'post_tt', '', ...
    'pre_firing_rate', 0, ...
    'post_firing_rate', 0, ...
    'connection_type_id', 0, ...
    'connection_type', '', ...
    'auto_type_id', 0, ...
    'auto_type', '', ...
    'initial_decision_id', 0, ...
    'initial_decision', '', ...
    'switched_from_auto', false, ...
    'same_tt_pair', false, ...
    'rebound_override', false, ...
    'verified', false, ...
    'last_reviewed', '', ...
    'presynaptic_cell_label', 0 ...
    );
end

function write_session_connection_exports(curated_session, curated_file)

accepted_connections = curated_session.accepted_connections;
if isempty(accepted_connections)
    accepted_connections = repmat(accepted_connection_template(), 0, 1);
end

accepted_table = struct2table(accepted_connections);
csv_file = strrep(curated_file, '.mat', '-AcceptedConnections.csv');
summary_file = strrep(curated_file, '.mat', '-Summary.txt');

writetable(accepted_table, csv_file);

fid = fopen(summary_file, 'w');
fprintf(fid, 'Session: %s\n', curated_session.session_label);
fprintf(fid, 'Curated on: %s\n', curated_session.curated_on);
fprintf(fid, 'Pairs evaluated: %d\n', curated_session.total_pairs_evaluated);
fprintf(fid, 'Identified candidate connections: %d\n', curated_session.total_connections_identified);
fprintf(fid, 'Verified candidate connections: %d\n', curated_session.total_connections_verified);
fprintf(fid, 'Accepted connections: %d\n', curated_session.accepted_connection_count);
fprintf(fid, 'Accepted plot count: %d\n', curated_session.accepted_plot_count);
fprintf(fid, 'Accepted plot folder: %s\n', empty_as_none(curated_session.accepted_plot_folder));
fprintf(fid, 'Accepted excitatory cross-TT pairs: %d\n', size(curated_session.connectivity.accepted_exc_pairs_crossTT, 1));
fprintf(fid, 'Accepted inhibitory cross-TT pairs: %d\n', size(curated_session.connectivity.accepted_inh_pairs_crossTT, 1));
fprintf(fid, 'Cell labels: excitatory %d | inhibitory %d | mixed %d | unlabeled %d\n', ...
    sum(curated_session.final_cell_classification == 1), ...
    sum(curated_session.final_cell_classification == 2), ...
    sum(curated_session.final_cell_classification == 3), ...
    sum(curated_session.final_cell_classification == 0));
fclose(fid);
end

function batch_summary = build_batch_summary(review_files, process_mode, skipped_empty_summary)

if nargin < 3 || isempty(skipped_empty_summary)
    skipped_empty_summary = struct('session_count', 0, 'pairs_evaluated', 0);
end

session_summaries = repmat(struct( ...
    'session_index', 0, ...
    'session_label', '', ...
    'review_file', '', ...
    'pairs_evaluated', 0, ...
    'identified_connections', 0, ...
    'verified_connections', 0, ...
    'accepted_excitatory_pairs', 0, ...
    'accepted_inhibitory_pairs', 0, ...
    'excitatory_cells', 0, ...
    'inhibitory_cells', 0, ...
    'mixed_cells', 0, ...
    'unlabeled_cells', 0, ...
    'review_complete', false, ...
    'writeback_complete', false ...
    ), numel(review_files), 1);

for idx = 1:numel(review_files)
    [review_session, ~] = load_review_session_file(review_files{idx});
    review_session = synchronize_review_session_status(review_session);
    cell_classification = compute_reviewed_cell_labels(review_session);

    session_summaries(idx).session_index = review_session.session_index;
    session_summaries(idx).session_label = review_session.session_label;
    session_summaries(idx).review_file = review_files{idx};
    session_summaries(idx).pairs_evaluated = review_session.total_pairs_evaluated;
    session_summaries(idx).identified_connections = review_session.total_connections_identified;
    session_summaries(idx).verified_connections = review_session.total_connections_verified;
    session_summaries(idx).accepted_excitatory_pairs = sum([review_session.review_pairs.current_decision] == 1);
    session_summaries(idx).accepted_inhibitory_pairs = sum([review_session.review_pairs.current_decision] == 2);
    session_summaries(idx).excitatory_cells = sum(cell_classification == 1);
    session_summaries(idx).inhibitory_cells = sum(cell_classification == 2);
    session_summaries(idx).mixed_cells = sum(cell_classification == 3);
    session_summaries(idx).unlabeled_cells = sum(cell_classification == 0);
    session_summaries(idx).review_complete = review_session.review_complete;
    session_summaries(idx).writeback_complete = review_session.writeback_complete;
end

batch_summary = struct();
batch_summary.generated_on = char(datetime('now'));
batch_summary.process_mode = process_mode;
batch_summary.review_file_count = numel(review_files);
batch_summary.skipped_empty_session_count = skipped_empty_summary.session_count;
batch_summary.skipped_empty_pairs_evaluated = skipped_empty_summary.pairs_evaluated;
batch_summary.total_pairs_evaluated = sum([session_summaries.pairs_evaluated]);
batch_summary.total_connections_tested = batch_summary.total_pairs_evaluated;
batch_summary.total_pairs_evaluated_including_empty = ...
    batch_summary.total_pairs_evaluated + batch_summary.skipped_empty_pairs_evaluated;
batch_summary.total_connections_identified = sum([session_summaries.identified_connections]);
batch_summary.total_connections_verified = sum([session_summaries.verified_connections]);
batch_summary.total_accepted_excitatory_pairs = sum([session_summaries.accepted_excitatory_pairs]);
batch_summary.total_accepted_inhibitory_pairs = sum([session_summaries.accepted_inhibitory_pairs]);
batch_summary.total_excitatory_cells = sum([session_summaries.excitatory_cells]);
batch_summary.total_inhibitory_cells = sum([session_summaries.inhibitory_cells]);
batch_summary.total_mixed_cells = sum([session_summaries.mixed_cells]);
batch_summary.total_unlabeled_cells = sum([session_summaries.unlabeled_cells]);
batch_summary.all_sessions_review_complete = all([session_summaries.review_complete]);
batch_summary.all_sessions_written_back = all([session_summaries.writeback_complete]);
batch_summary.session_summaries = session_summaries;
end

function latest_summary_file = save_batch_summary(batch_summary, summary_root)

ensure_dir(summary_root);

latest_summary_file = fullfile(summary_root, 'CCG_ReviewSummary_Latest.txt');
latest_summary_mat = fullfile(summary_root, 'CCG_ReviewSummary_Latest.mat');

save(latest_summary_mat, 'batch_summary', '-v7.3');

fid = fopen(latest_summary_file, 'w');
fprintf(fid, 'Generated: %s\n', batch_summary.generated_on);
fprintf(fid, 'Mode: %s\n', process_mode_label(batch_summary.process_mode));
fprintf(fid, 'Review files: %d\n', batch_summary.review_file_count);
fprintf(fid, 'Total pairs evaluated: %d\n', batch_summary.total_pairs_evaluated);
fprintf(fid, 'Skipped empty sessions: %d\n', batch_summary.skipped_empty_session_count);
fprintf(fid, 'Pairs evaluated in skipped empty sessions: %d\n', batch_summary.skipped_empty_pairs_evaluated);
fprintf(fid, 'Total pairs evaluated including skipped empty sessions: %d\n', ...
    batch_summary.total_pairs_evaluated_including_empty);
fprintf(fid, 'Total identified candidate connections: %d\n', batch_summary.total_connections_identified);
fprintf(fid, 'Total verified candidate connections: %d\n', batch_summary.total_connections_verified);
fprintf(fid, 'Accepted excitatory pairs: %d\n', batch_summary.total_accepted_excitatory_pairs);
fprintf(fid, 'Accepted inhibitory pairs: %d\n', batch_summary.total_accepted_inhibitory_pairs);
fprintf(fid, 'Cell labels: excitatory %d | inhibitory %d | mixed %d | unlabeled %d\n', ...
    batch_summary.total_excitatory_cells, batch_summary.total_inhibitory_cells, ...
    batch_summary.total_mixed_cells, batch_summary.total_unlabeled_cells);
fprintf(fid, 'All sessions review-complete: %s\n', ternary(batch_summary.all_sessions_review_complete, 'yes', 'no'));
fprintf(fid, 'All sessions written back: %s\n\n', ternary(batch_summary.all_sessions_written_back, 'yes', 'no'));
fprintf(fid, 'Per-session summary:\n');

for idx = 1:numel(batch_summary.session_summaries)
    entry = batch_summary.session_summaries(idx);
    fprintf(fid, ['- %s | evaluated %d | identified %d | verified %d | accepted exc %d | accepted inh %d | ' ...
        'cells exc %d inh %d mixed %d unlabeled %d | writeback %s\n'], ...
        entry.session_label, entry.pairs_evaluated, entry.identified_connections, entry.verified_connections, ...
        entry.accepted_excitatory_pairs, entry.accepted_inhibitory_pairs, ...
        entry.excitatory_cells, entry.inhibitory_cells, entry.mixed_cells, ...
        entry.unlabeled_cells, ternary(entry.writeback_complete, 'yes', 'no'));
end
fclose(fid);

if batch_summary.all_sessions_written_back
    stamp = datestr(now, 'yyyymmdd_HHMMSS');
    save(fullfile(summary_root, sprintf('CCG_ReviewSummary_%s.mat', stamp)), 'batch_summary', '-v7.3');
end
end

function global_summary = build_global_curated_summary(all_cells_file, curated_session_root)

loaded_all = load(all_cells_file, 'All_Cells_combined');
if ~isfield(loaded_all, 'All_Cells_combined')
    error('All_Cells_combined was not found in %s.', all_cells_file);
end

All_Cells_combined = loaded_all.All_Cells_combined;
n_sessions = numel(All_Cells_combined);

session_summaries = repmat(struct( ...
    'session_index', 0, ...
    'session_label', '', ...
    'animal', [], ...
    'day', [], ...
    'unit_count', 0, ...
    'identified_connections', 0, ...
    'verified_connections', 0, ...
    'accepted_excitatory_pairs', 0, ...
    'accepted_inhibitory_pairs', 0, ...
    'accepted_connection_count', 0, ...
    'excitatory_cells', 0, ...
    'inhibitory_cells', 0, ...
    'mixed_cells', 0, ...
    'unlabeled_cells', 0, ...
    'has_curated_file', false, ...
    'zero_only_session', false, ...
    'curated_file', '' ...
    ), n_sessions, 1);

classifications = cell(n_sessions, 1);

for idx = 1:n_sessions
    session_entry = All_Cells_combined(idx);
    unit_count = infer_all_cells_session_unit_count_for_summary(session_entry);
    classifications{idx} = zeros(unit_count, 1);
    session_summaries(idx).session_index = idx;
    session_summaries(idx).session_label = build_global_session_label(session_entry, idx);
    session_summaries(idx).animal = get_struct_field_or_empty(session_entry, 'animal');
    session_summaries(idx).day = get_struct_field_or_empty(session_entry, 'day');
    session_summaries(idx).unit_count = unit_count;
end

listing = dir(fullfile(curated_session_root, '*-Curated.mat'));
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

    base_count = numel(classifications{session_index});
    curated_count = infer_curated_session_unit_count_for_summary(curated_session);
    final_count = max([base_count, curated_count, numel(curated_session.final_cell_classification)]);
    classification = zeros(final_count, 1);
    source_classification = curated_session.final_cell_classification(:);

    if ~isempty(source_classification)
        classification(1:numel(source_classification)) = source_classification;
    end

    classifications{session_index} = classification;
    session_summaries(session_index).unit_count = final_count;
    session_summaries(session_index).session_label = curated_session.session_label;
    session_summaries(session_index).animal = curated_session.animal;
    session_summaries(session_index).day = curated_session.day;
    session_summaries(session_index).identified_connections = get_struct_field_or_default(curated_session, 'total_connections_identified', 0);
    session_summaries(session_index).verified_connections = get_struct_field_or_default(curated_session, 'total_connections_verified', 0);
    session_summaries(session_index).accepted_excitatory_pairs = count_curated_connections(curated_session, 1);
    session_summaries(session_index).accepted_inhibitory_pairs = count_curated_connections(curated_session, 2);
    session_summaries(session_index).accepted_connection_count = numel(get_struct_field_or_default(curated_session, 'accepted_connections', []));
    session_summaries(session_index).has_curated_file = true;
    session_summaries(session_index).curated_file = fullfile(listing(idx).folder, listing(idx).name);
end

for idx = 1:n_sessions
    classification = classifications{idx};
    session_summaries(idx).excitatory_cells = sum(classification == 1);
    session_summaries(idx).inhibitory_cells = sum(classification == 2);
    session_summaries(idx).mixed_cells = sum(classification == 3);
    session_summaries(idx).unlabeled_cells = sum(classification == 0);
    session_summaries(idx).zero_only_session = all(classification == 0);
end

global_summary = struct();
global_summary.generated_on = char(datetime('now'));
global_summary.all_cells_file = all_cells_file;
global_summary.curated_session_root = curated_session_root;
global_summary.total_sessions = n_sessions;
global_summary.curated_session_count = sum([session_summaries.has_curated_file]);
global_summary.sessions_without_curated_file = n_sessions - global_summary.curated_session_count;
global_summary.zero_only_session_count = sum([session_summaries.zero_only_session]);
global_summary.total_cells = sum([session_summaries.unit_count]);
global_summary.total_identified_connections = sum([session_summaries.identified_connections]);
global_summary.total_verified_connections = sum([session_summaries.verified_connections]);
global_summary.total_accepted_excitatory_pairs = sum([session_summaries.accepted_excitatory_pairs]);
global_summary.total_accepted_inhibitory_pairs = sum([session_summaries.accepted_inhibitory_pairs]);
global_summary.total_accepted_connections = sum([session_summaries.accepted_connection_count]);
global_summary.total_excitatory_cells = sum([session_summaries.excitatory_cells]);
global_summary.total_inhibitory_cells = sum([session_summaries.inhibitory_cells]);
global_summary.total_mixed_cells = sum([session_summaries.mixed_cells]);
global_summary.total_unlabeled_cells = sum([session_summaries.unlabeled_cells]);
global_summary.skipped_curated_entries = skipped_sessions(:)';
global_summary.session_summaries = session_summaries;
end

function latest_summary_file = save_global_curated_summary(global_summary, summary_root)

ensure_dir(summary_root);

latest_summary_file = fullfile(summary_root, 'CCG_ReviewSummary_Latest.txt');
latest_summary_mat = fullfile(summary_root, 'CCG_ReviewSummary_Latest.mat');
global_summary_file = fullfile(summary_root, 'CCG_GlobalCuratedSummary_Latest.txt');
global_summary_mat = fullfile(summary_root, 'CCG_GlobalCuratedSummary_Latest.mat');
global_session_csv = fullfile(summary_root, 'CCG_GlobalSessionSummary.csv');

save(latest_summary_mat, 'global_summary', '-v7.3');
save(global_summary_mat, 'global_summary', '-v7.3');
writetable(struct2table(global_summary.session_summaries), global_session_csv);

write_global_summary_text(latest_summary_file, global_summary);
write_global_summary_text(global_summary_file, global_summary);
end

function write_global_summary_text(summary_file, global_summary)

fid = fopen(summary_file, 'w');
fprintf(fid, 'Generated: %s\n', global_summary.generated_on);
fprintf(fid, 'Summary type: Global curated summary\n');
fprintf(fid, 'All_Cells file: %s\n', global_summary.all_cells_file);
fprintf(fid, 'Curated session folder: %s\n', global_summary.curated_session_root);
fprintf(fid, 'Total sessions in All_Cells_combined: %d\n', global_summary.total_sessions);
fprintf(fid, 'Sessions with curated files: %d\n', global_summary.curated_session_count);
fprintf(fid, 'Sessions without curated files: %d\n', global_summary.sessions_without_curated_file);
fprintf(fid, 'Zero-only sessions: %d\n', global_summary.zero_only_session_count);
fprintf(fid, 'Total cells across all sessions: %d\n', global_summary.total_cells);
fprintf(fid, 'Total identified candidate connections across curated sessions: %d\n', global_summary.total_identified_connections);
fprintf(fid, 'Total verified candidate connections across curated sessions: %d\n', global_summary.total_verified_connections);
fprintf(fid, 'Accepted excitatory pairs across curated sessions: %d\n', global_summary.total_accepted_excitatory_pairs);
fprintf(fid, 'Accepted inhibitory pairs across curated sessions: %d\n', global_summary.total_accepted_inhibitory_pairs);
fprintf(fid, 'Accepted connections total: %d\n', global_summary.total_accepted_connections);
fprintf(fid, 'Cell labels: excitatory %d | inhibitory %d | mixed %d | unlabeled %d\n', ...
    global_summary.total_excitatory_cells, global_summary.total_inhibitory_cells, ...
    global_summary.total_mixed_cells, global_summary.total_unlabeled_cells);

if ~isempty(global_summary.skipped_curated_entries)
    fprintf(fid, 'Skipped curated entries: %s\n', mat2str(global_summary.skipped_curated_entries));
end

fprintf(fid, '\nPer-session summary:\n');

for idx = 1:numel(global_summary.session_summaries)
    entry = global_summary.session_summaries(idx);
    fprintf(fid, ['- %s | cells total %d | identified %d | verified %d | accepted exc %d | accepted inh %d | ' ...
        'cells exc %d inh %d mixed %d unlabeled %d | curated %s\n'], ...
        entry.session_label, entry.unit_count, entry.identified_connections, entry.verified_connections, ...
        entry.accepted_excitatory_pairs, entry.accepted_inhibitory_pairs, ...
        entry.excitatory_cells, entry.inhibitory_cells, entry.mixed_cells, ...
        entry.unlabeled_cells, ternary(entry.has_curated_file, 'yes', 'no'));
end

fclose(fid);
end

function unit_count = infer_all_cells_session_unit_count_for_summary(session_entry)

unit_count = 0;

if isfield(session_entry, 'classific_firingRate') && ~isempty(session_entry.classific_firingRate)
    unit_count = numel(session_entry.classific_firingRate);
elseif isfield(session_entry, 'CCGbased_CellClassfication') && ~isempty(session_entry.CCGbased_CellClassfication)
    unit_count = numel(session_entry.CCGbased_CellClassfication);
end
end

function unit_count = infer_curated_session_unit_count_for_summary(curated_session)

unit_count = 0;

if isfield(curated_session, 'n_units') && ~isempty(curated_session.n_units)
    unit_count = double(curated_session.n_units);
elseif isfield(curated_session, 'firing_rates') && ~isempty(curated_session.firing_rates)
    unit_count = numel(curated_session.firing_rates);
elseif isfield(curated_session, 'final_cell_classification') && ~isempty(curated_session.final_cell_classification)
    unit_count = numel(curated_session.final_cell_classification);
end
end

function session_label = build_global_session_label(session_entry, session_index)

animal = get_struct_field_or_empty(session_entry, 'animal');
day = get_struct_field_or_empty(session_entry, 'day');
session_label = sprintf('Animal %s Day %s Session %d', ...
    stringify_animal_value(animal), stringify_animal_value(day), session_index);
end

function value = get_struct_field_or_empty(in_struct, field_name)

if isfield(in_struct, field_name)
    value = in_struct.(field_name);
else
    value = [];
end
end

function value = get_struct_field_or_default(in_struct, field_name, default_value)

if isfield(in_struct, field_name) && ~isempty(in_struct.(field_name))
    value = in_struct.(field_name);
else
    value = default_value;
end
end

function count = count_curated_connections(curated_session, connection_type_id)

count = 0;

if ~isfield(curated_session, 'accepted_connections') || isempty(curated_session.accepted_connections)
    return
end

count = sum([curated_session.accepted_connections.connection_type_id] == connection_type_id);
end

function update_global_curated_connectivity(curated_session_root, curated_root)

ensure_dir(curated_root);

listing = dir(fullfile(curated_session_root, '*-Curated.mat'));
if isempty(listing)
    return
end

all_connections = repmat(accepted_connection_template(), 0, 1);
session_summary = repmat(struct( ...
    'session_index', 0, ...
    'session_label', '', ...
    'animal', 0, ...
    'day', 0, ...
    'accepted_connection_count', 0, ...
    'total_connections_tested', 0, ...
    'total_connections_identified', 0, ...
    'total_connections_verified', 0, ...
    'excitatory_cells', 0, ...
    'inhibitory_cells', 0, ...
    'mixed_cells', 0, ...
    'unlabeled_cells', 0, ...
    'curated_file', '' ...
    ), numel(listing), 1);

for idx = 1:numel(listing)
    loaded = load(fullfile(listing(idx).folder, listing(idx).name), 'curated_session');
    curated_session = loaded.curated_session;

    if ~isempty(curated_session.accepted_connections)
        all_connections = [all_connections; curated_session.accepted_connections]; %#ok<AGROW>
    end

    session_summary(idx).session_index = curated_session.session_index;
    session_summary(idx).session_label = curated_session.session_label;
    session_summary(idx).animal = curated_session.animal;
    session_summary(idx).day = curated_session.day;
    session_summary(idx).accepted_connection_count = curated_session.accepted_connection_count;
    session_summary(idx).total_connections_tested = curated_session.total_connections_tested;
    session_summary(idx).total_connections_identified = curated_session.total_connections_identified;
    session_summary(idx).total_connections_verified = curated_session.total_connections_verified;
    session_summary(idx).excitatory_cells = sum(curated_session.final_cell_classification == 1);
    session_summary(idx).inhibitory_cells = sum(curated_session.final_cell_classification == 2);
    session_summary(idx).mixed_cells = sum(curated_session.final_cell_classification == 3);
    session_summary(idx).unlabeled_cells = sum(curated_session.final_cell_classification == 0);
    session_summary(idx).curated_file = fullfile(listing(idx).folder, listing(idx).name);
end

save(fullfile(curated_root, 'CCG_AcceptedConnections_AllSessions.mat'), ...
    'all_connections', 'session_summary', '-v7.3');

writetable(struct2table(all_connections), fullfile(curated_root, 'CCG_AcceptedConnections_AllSessions.csv'));
writetable(struct2table(session_summary), fullfile(curated_root, 'CCG_CuratedSessionSummary.csv'));
end

function label = type_label(type_id)

switch type_id
    case 1
        label = 'excitatory';
    case 2
        label = 'inhibitory';
    case 3
        label = 'both';
    otherwise
        label = 'none';
end
end

function label = decision_label(type_id)

switch type_id
    case 1
        label = 'accepted excitatory';
    case 2
        label = 'accepted inhibitory';
    otherwise
        label = 'discarded';
end
end

function plot_type = effective_plot_type(pair)

if pair.current_decision == 1 || pair.current_decision == 2
    plot_type = pair.current_decision;
elseif pair.initial_decision == 1 || pair.initial_decision == 2
    plot_type = pair.initial_decision;
elseif pair.rebound_override
    plot_type = 2;
elseif pair.auto_type == 2
    plot_type = 2;
else
    plot_type = 1;
end
end

function out = ternary(condition, true_value, false_value)

if condition
    out = true_value;
else
    out = false_value;
end
end

function out = on_off(condition)

if condition
    out = 'on';
else
    out = 'off';
end
end

function txt = empty_as_none(txt)

if isempty(txt)
    txt = 'none';
end
end

function label = process_mode_label(process_mode)

switch process_mode
    case 'export_plots_only'
        label = 'Export accepted plots only';
    case 'reprocess_all'
        label = 'Reprocess all';
    otherwise
        label = 'Only new/unverified';
end
end

function ensure_dir(folder_path)

if ~exist(folder_path, 'dir')
    mkdir(folder_path)
end
end

function signature = make_candidate_signature(review_pairs)

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

function fig_handle = render_presynaptic_cell_view(review_session, pre_idx, fig_handle)

review_session = synchronize_review_session_status(review_session);
tt_label = review_session.tt_files{pre_idx};

[acg_counts, acg_lags_s] = CrossCorrel(review_session.spike_times{pre_idx}, review_session.spike_times{pre_idx}, ...
    0.001, [-0.100 0.100]);
acg_lags_ms = acg_lags_s(:) * 1000;

map_entries = resolve_presynaptic_map_entries(review_session, pre_idx);
nPanels = 1 + numel(map_entries);
nCols = min(3, max(2, ceil(sqrt(nPanels))));
nRows = ceil(nPanels / nCols);

if isempty(fig_handle) || ~ishandle(fig_handle)
    fig_handle = figure( ...
        'Name', sprintf('Presynaptic Cell - %s', tt_label), ...
        'Units', 'normalized', ...
        'Position', [0.08 0.08 0.84 0.82], ...
        'Color', 'w', ...
        'NumberTitle', 'off' ...
        );
else
    figure(fig_handle);
    clf(fig_handle);
    set(fig_handle, 'Name', sprintf('Presynaptic Cell - %s', tt_label));
end

subplot(nRows, nCols, 1, 'Parent', fig_handle);
bar(acg_lags_ms, double(acg_counts(:)), 1, 'FaceColor', [0.3 0.3 0.3], 'EdgeColor', 'none');
xlim([-100 100]);
xlabel('time [ms]');
ylabel('Count');
title(sprintf('Autocorrelogram | %s | FR %.3f Hz', ...
    tt_label, review_session.firing_rates(pre_idx)));

for idx = 1:numel(map_entries)
    subplot(nRows, nCols, idx + 1, 'Parent', fig_handle);
    entry = map_entries(idx);

    if strcmp(entry.display_mode, 'map') && ~isempty(entry.map_matrix)
        finite_mask = isfinite(entry.map_matrix);
        imagesc(entry.map_matrix, 'AlphaData', double(finite_mask));
        set(gca, 'Color', [0.82 0.82 0.82], 'YDir', 'normal');
        axis image off;
        colormap(gca, jet(256));
        if any(finite_mask(:))
            colorbar;
        else
            text(0.5, 0.5, 'No mapped spikes', 'Units', 'normalized', ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        end
        if isempty(entry.status_text)
            title(sprintf('%s | %s', entry.session_dir, entry.source_label), 'Interpreter', 'none');
        else
            title(sprintf('%s | %s\n%s', entry.session_dir, entry.source_label, entry.status_text), ...
                'Interpreter', 'none');
        end
    elseif ~isempty(entry.image_path) && exist(entry.image_path, 'file') == 2
        image_data = imread(entry.image_path);
        image(image_data);
        axis image off;
        title(sprintf('%s | %s', entry.session_dir, entry.source_label), 'Interpreter', 'none');
    else
        axis off;
        text(0.5, 0.58, entry.session_dir, 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Interpreter', 'none');
        text(0.5, 0.42, entry.status_text, 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'Interpreter', 'none');
        title('Spatial Map');
    end
end
end

function map_entries = resolve_presynaptic_map_entries(review_session, pre_idx)

behavior_dirs = get_behavior_session_dirs(review_session);
tt_file = review_session.tt_files{pre_idx};
[~, tt_stem, tt_ext] = fileparts(tt_file);
if isempty(tt_stem)
    tt_stem = tt_file;
else
    tt_stem = [tt_stem, tt_ext];
    [~, tt_stem] = fileparts(tt_stem);
end

if isempty(behavior_dirs)
    map_entries = struct('session_dir', 'Behavior Sessions', 'image_path', '', ...
        'source_label', '', 'status_text', 'No behavioral subsessions were available.', ...
        'display_mode', 'text', 'map_matrix', []);
    return
end

map_entries = repmat(struct( ...
    'session_dir', '', ...
    'image_path', '', ...
    'source_label', '', ...
    'status_text', '', ...
    'display_mode', 'text', ...
    'map_matrix', [] ...
    ), numel(behavior_dirs), 1);

if has_exact_behavior_spike_data(review_session, behavior_dirs)
    map_entries = build_exact_presynaptic_map_entries(review_session, pre_idx, behavior_dirs);
    return
end

missing_paths = {};
missing_idx = [];

for idx = 1:numel(behavior_dirs)
    session_dir = behavior_dirs{idx};
    session_path = ensure_trailing_filesep(fullfile(review_session.mainDir, session_dir));
    [image_path, source_label] = find_existing_presynaptic_map_image(session_path, tt_stem);

    map_entries(idx).session_dir = session_dir;
    map_entries(idx).image_path = image_path;
    map_entries(idx).source_label = source_label;

    if isempty(image_path)
        map_entries(idx).status_text = 'Map not cached yet.';
        map_entries(idx).display_mode = 'text';
        missing_paths{end+1,1} = session_path; %#ok<AGROW>
        missing_idx(end+1,1) = idx; %#ok<AGROW>
    else
        map_entries(idx).status_text = '';
        map_entries(idx).display_mode = 'image';
    end
end

if isempty(missing_paths)
    return
end

compute_error = '';
try
    compute_presynaptic_maps(missing_paths, tt_file);
catch exception
    compute_error = exception.message;
end

for idx = 1:numel(missing_idx)
    entry_idx = missing_idx(idx);
    session_path = missing_paths{idx};
    [image_path, source_label] = find_existing_presynaptic_map_image(session_path, tt_stem);

    map_entries(entry_idx).image_path = image_path;
    map_entries(entry_idx).source_label = source_label;

    if isempty(image_path)
        if isempty(compute_error)
            map_entries(entry_idx).status_text = 'Map could not be generated.';
        else
            map_entries(entry_idx).status_text = compute_error;
        end
        map_entries(entry_idx).display_mode = 'text';
    else
        map_entries(entry_idx).status_text = '';
        map_entries(entry_idx).display_mode = 'image';
    end
end
end

function tf = has_exact_behavior_spike_data(review_session, behavior_dirs)

tf = isfield(review_session, 'behavior_spike_times') && ...
    ~isempty(review_session.behavior_spike_times) && ...
    numel(review_session.behavior_spike_times) >= numel(behavior_dirs);
end

function map_entries = build_exact_presynaptic_map_entries(review_session, pre_idx, behavior_dirs)

map_entries = repmat(struct( ...
    'session_dir', '', ...
    'image_path', '', ...
    'source_label', '', ...
    'status_text', '', ...
    'display_mode', 'text', ...
    'map_matrix', [] ...
    ), numel(behavior_dirs), 1);

for idx = 1:numel(behavior_dirs)
    session_dir = behavior_dirs{idx};
    session_path = ensure_trailing_filesep(fullfile(review_session.mainDir, session_dir));

    map_entries(idx).session_dir = session_dir;
    map_entries(idx).source_label = 'exact review spikes';

    if idx > numel(review_session.behavior_spike_times) || isempty(review_session.behavior_spike_times{idx})
        map_entries(idx).status_text = 'No per-session spike data available.';
        continue
    end

    spike_times_this_session = review_session.behavior_spike_times{idx};
    if pre_idx > numel(spike_times_this_session)
        map_entries(idx).status_text = 'Presynaptic cell was not present in this subsession.';
        continue
    end

    try
        [rate_map, map_status] = compute_exact_presynaptic_map(session_path, spike_times_this_session{pre_idx});
    catch exception
        rate_map = [];
        map_status = exception.message;
    end

    if isempty(rate_map)
        map_entries(idx).status_text = map_status;
        continue
    end

    map_entries(idx).map_matrix = rate_map;
    map_entries(idx).status_text = map_status;
    map_entries(idx).display_mode = 'map';
end
end

function [rate_map, status_text] = compute_exact_presynaptic_map(session_path, spike_times_raw)

ensure_cluster_cutting_paths();

rate_map = [];
status_text = '';

spike_times_raw = double(spike_times_raw(:));
spike_times_raw = spike_times_raw(isfinite(spike_times_raw));
if isempty(spike_times_raw)
    status_text = 'No spike timestamps available.';
    return
end

tracking_file = local_resolve_review_tracking_file(session_path);
if isempty(tracking_file)
    status_text = 'Tracking file not found.';
    return
end

[t, x, y] = getTrackingData(tracking_file, 7, 1);
t = double(t(:));
x = double(x(:));
y = double(y(:));

finite_mask = isfinite(t) & isfinite(x) & isfinite(y);
t = t(finite_mask);
x = x(finite_mask);
y = y(finite_mask);

if isempty(t)
    status_text = 'Tracking data is empty.';
    return
end

[t, sort_idx] = sort(t);
x = x(sort_idx);
y = y(sort_idx);

if numel(t) < 2
    status_text = 'Tracking data is too short.';
    return
end

[x, y] = local_auto_center_and_scale_review_tracking(x, y, 50, 50);

alignment_info = local_infer_review_spike_time_transform(spike_times_raw, t, 1);
spike_times_sec = local_apply_review_spike_time_transform(spike_times_raw, alignment_info);

in_range_mask = spike_times_sec >= t(1) & spike_times_sec < t(end);
spike_times_in_range = spike_times_sec(in_range_mask);

if isempty(spike_times_in_range)
    status_text = sprintf('No spikes in tracking window | factor %.6g | offset %.3f s | overlap %.2f', ...
        alignment_info.factor, alignment_info.offsetSec, alignment_info.overlapFraction);
    return
end

[spkx, spky] = spikePos(spike_times_in_range, x, y, t, 'interpolate', 0.5);
assigned_mask = isfinite(spkx) & isfinite(spky);
assigned_count = sum(assigned_mask);
total_count = numel(spike_times_raw);
in_range_count = numel(spike_times_in_range);

if assigned_count == 0
    status_text = sprintf('No spikes mapped (0/%d assigned, %d in range) | factor %.6g | offset %.3f s | overlap %.2f', ...
        total_count, in_range_count, alignment_info.factor, alignment_info.offsetSec, alignment_info.overlapFraction);
    return
end

spike_times_in_range = spike_times_in_range(assigned_mask);
spkx = spkx(assigned_mask);
spky = spky(assigned_mask);

geometry = local_build_review_map_geometry(50, 50, 50);
time_map = local_find_review_time_map(x, y, t, geometry.xLimits, geometry.yLimits, ...
    geometry.binWidthX, geometry.binWidthY, geometry.binsX, geometry.binsY);

if ~any(time_map(:) > 0)
    status_text = 'No visited bins in tracking data.';
    return
end

rate_map = local_review_rate_map(spike_times_in_range, spkx, spky, x, y, t, 5, geometry.xAxis, geometry.yAxis);
rate_map(time_map == 0) = NaN;
rate_map(~isfinite(rate_map)) = NaN;

status_text = sprintf('Assigned %d/%d spikes (%d in range) | factor %.6g | offset %.3f s | overlap %.2f', ...
    assigned_count, total_count, in_range_count, alignment_info.factor, alignment_info.offsetSec, alignment_info.overlapFraction);
end

function tracking_file = local_resolve_review_tracking_file(session_path)

tracking_file = '';
candidate_names = {'vt1.nvt', 'VT1.nvt', 'VT.nvt', 'vt.nvt'};

for idx = 1:numel(candidate_names)
    candidate_path = fullfile(session_path, candidate_names{idx});
    if exist(candidate_path, 'file') == 2
        tracking_file = candidate_path;
        return
    end
end

listing = dir(fullfile(session_path, '*.nvt'));
if ~isempty(listing)
    tracking_file = fullfile(session_path, listing(1).name);
end
end

function alignment_info = local_infer_review_spike_time_transform(ts_raw, tracking_times_sec, default_factor)

if nargin < 3 || isempty(default_factor)
    default_factor = 1;
end

ts_raw = double(ts_raw(:));
tracking_times_sec = double(tracking_times_sec(:));
tracking_times_sec = tracking_times_sec(isfinite(tracking_times_sec));

alignment_info = struct('factor', default_factor, 'offsetSec', 0, ...
    'score', -inf, 'overlapFraction', 0, 'rangeError', inf, 'centerError', inf);

if isempty(ts_raw) || isempty(tracking_times_sec)
    return
end

t_min = min(tracking_times_sec);
t_max = max(tracking_times_sec);
t_range = max(t_max - t_min, eps);
t_center = (t_min + t_max) / 2;
candidate_factors = unique([default_factor, 1, 1e-4, 1e-6, 1e-3, 1e-8], 'stable');

for idx_factor = 1:numel(candidate_factors)
    factor = candidate_factors(idx_factor);
    ts_scaled = ts_raw * factor;
    valid_mask = isfinite(ts_scaled);
    if ~any(valid_mask)
        continue
    end

    ts_valid = ts_scaled(valid_mask);
    ts_min = min(ts_valid);
    ts_max = max(ts_valid);
    ts_center = (ts_min + ts_max) / 2;
    candidate_offsets = unique([0, t_min - ts_min, t_max - ts_max, t_center - ts_center]);

    for idx_offset = 1:numel(candidate_offsets)
        offset_sec = candidate_offsets(idx_offset);
        ts_candidate = ts_scaled + offset_sec;
        [score, overlap_fraction, range_error, center_error] = ...
            local_score_review_timestamp_alignment(ts_candidate, t_min, t_max, t_range, t_center);

        if score > alignment_info.score || ...
                (abs(score - alignment_info.score) < 1e-12 && overlap_fraction > alignment_info.overlapFraction) || ...
                (abs(score - alignment_info.score) < 1e-12 && abs(overlap_fraction - alignment_info.overlapFraction) < 1e-12 && range_error < alignment_info.rangeError) || ...
                (abs(score - alignment_info.score) < 1e-12 && abs(range_error - alignment_info.rangeError) < 1e-12 && center_error < alignment_info.centerError)
            alignment_info.factor = factor;
            alignment_info.offsetSec = offset_sec;
            alignment_info.score = score;
            alignment_info.overlapFraction = overlap_fraction;
            alignment_info.rangeError = range_error;
            alignment_info.centerError = center_error;
        end
    end
end
end

function [score, overlap_fraction, range_error, center_error] = ...
        local_score_review_timestamp_alignment(ts_seconds, t_min, t_max, t_range, t_center)

valid_mask = isfinite(ts_seconds);
if ~any(valid_mask)
    overlap_fraction = 0;
    range_error = inf;
    center_error = inf;
    score = -inf;
    return
end

ts_valid = ts_seconds(valid_mask);
range_margin = max(1, 0.05 * t_range);
in_range_mask = ts_valid >= (t_min - range_margin) & ts_valid <= (t_max + range_margin);
overlap_fraction = mean(in_range_mask);

ts_min = min(ts_valid);
ts_max = max(ts_valid);
ts_center = (ts_min + ts_max) / 2;
range_error = abs((ts_max - ts_min) - t_range) / t_range;
center_error = abs(ts_center - t_center) / t_range;
score = overlap_fraction - 0.05 * range_error - 0.02 * center_error;
end

function ts_seconds = local_apply_review_spike_time_transform(ts_raw, alignment_info)

ts_seconds = double(ts_raw(:)) * alignment_info.factor + alignment_info.offsetSec;
end

function geometry = local_build_review_map_geometry(track_width_cm, track_height_cm, bins_x)

bin_width_x = track_width_cm / bins_x;
bins_y = max(1, round(track_height_cm / bin_width_x));
bin_width_y = track_height_cm / bins_y;

x_limits = [-track_width_cm/2, track_width_cm/2];
y_limits = [-track_height_cm/2, track_height_cm/2];
x_axis = linspace(x_limits(1) + bin_width_x/2, x_limits(2) - bin_width_x/2, bins_x);
y_axis = linspace(y_limits(1) + bin_width_y/2, y_limits(2) - bin_width_y/2, bins_y);

geometry = struct('xAxis', x_axis, 'yAxis', y_axis, ...
    'xLimits', x_limits, 'yLimits', y_limits, ...
    'binWidthX', bin_width_x, 'binWidthY', bin_width_y, ...
    'binsX', bins_x, 'binsY', bins_y);
end

function [posx, posy] = local_auto_center_and_scale_review_tracking(posx, posy, track_width_cm, track_height_cm)

posx = posx(:);
posy = posy(:);

finite_mask = isfinite(posx) & isfinite(posy);
if ~any(finite_mask)
    error('Tracking data is empty after filtering.');
end

x_min = min(posx(finite_mask));
x_max = max(posx(finite_mask));
y_min = min(posy(finite_mask));
y_max = max(posy(finite_mask));

x_centre = (x_min + x_max) / 2;
y_centre = (y_min + y_max) / 2;

x_range = x_max - x_min;
y_range = y_max - y_min;
if x_range <= 0
    x_range = 1;
end
if y_range <= 0
    y_range = 1;
end

posx = (posx - x_centre) * (track_width_cm / x_range);
posy = (posy - y_centre) * (track_height_cm / y_range);
end

function time_map = local_find_review_time_map(posx, posy, post, x_limits, y_limits, ...
        bin_width_x, bin_width_y, bins_x, bins_y)

nan_inds = isnan(posx + posy + post);
posx(nan_inds) = [];
posy(nan_inds) = [];
post(nan_inds) = [];

duration = post(end) - post(1);
samp_dur = duration / length(posx);

time_map = zeros(bins_y, bins_x);
for ii = 1:bins_x
    x_low = x_limits(1) + (ii-1) * bin_width_x;
    x_high = x_low + bin_width_x;
    if ii == bins_x
        idx_x = find(posx >= x_low & posx <= x_high);
    else
        idx_x = find(posx >= x_low & posx < x_high);
    end

    for jj = 1:bins_y
        y_low = y_limits(1) + (jj-1) * bin_width_y;
        y_high = y_low + bin_width_y;
        if jj == bins_y
            idx_y = find(posy(idx_x) >= y_low & posy(idx_x) <= y_high);
        else
            idx_y = find(posy(idx_x) >= y_low & posy(idx_x) < y_high);
        end
        time_map(jj, ii) = length(idx_y);
    end
end

time_map = time_map * samp_dur;
end

function map = local_review_rate_map(ts, spkx, spky, posx, posy, post, h, x_axis, y_axis)

nan_inds = isnan(posx + posy + post);
posx(nan_inds) = [];
posy(nan_inds) = [];
post(nan_inds) = [];

invh = 1 / h;
map = zeros(length(y_axis), length(x_axis));
yy = 0;
for y = y_axis
    yy = yy + 1;
    xx = 0;
    for x = x_axis
        xx = xx + 1;
        map(yy, xx) = local_review_rate_estimator(ts, spkx, spky, x, y, invh, posx, posy, post);
    end
end
end

function rate_value = local_review_rate_estimator(ts, spkx, spky, x, y, invh, posx, posy, post)

conv_sum = sum(local_review_gaussian_kernel((spkx - x) * invh, (spky - y) * invh));
edge_corrector = trapz(post, local_review_gaussian_kernel((posx - x) * invh, (posy - y) * invh));
rate_value = (conv_sum / (edge_corrector + 0.1)) + 0.1;
end

function values = local_review_gaussian_kernel(x, y)

values = 0.15915494309190 * exp(-0.5 * (x .* x + y .* y));
end

function behavior_dirs = get_behavior_session_dirs(review_session)

behavior_dirs = {};

if isfield(review_session, 'behavior_session_dirs') && ~isempty(review_session.behavior_session_dirs)
    behavior_dirs = review_session.behavior_session_dirs(:);
    return
end

session_info_file = resolve_session_info_path();
if exist(session_info_file, 'file') ~= 2
    return
end

loaded = load(session_info_file, 'sessInfo');
if ~isfield(loaded, 'sessInfo')
    return
end

if review_session.session_index < 1 || review_session.session_index > numel(loaded.sessInfo)
    return
end

behavior_dirs = loaded.sessInfo(review_session.session_index).sessDirs(:);
end

function [image_path, source_label] = find_existing_presynaptic_map_image(session_path, tt_stem)

image_path = '';
source_label = '';

folder_listing = dir(fullfile(session_path, 'placeFieldImages*'));
folder_listing = folder_listing([folder_listing.isdir]);

if isempty(folder_listing)
    return
end

folder_names = {folder_listing.name};
folder_names = folder_names(~ismember(folder_names, {'.', '..'}));

is_review_cache = cellfun(@(name) contains(name, 'ccgReview'), folder_names);
[~, order] = sort(double(~is_review_cache));
folder_names = folder_names(order);

extensions = {'bmp', 'jpg', 'png'};
for idx_dir = 1:numel(folder_names)
    folder_name = folder_names{idx_dir};
    for idx_ext = 1:numel(extensions)
        candidate_path = fullfile(session_path, folder_name, [tt_stem, '.', extensions{idx_ext}]);
        if exist(candidate_path, 'file') == 2
            image_path = candidate_path;
            if contains(folder_name, 'ccgReview')
                source_label = 'review cache';
            else
                source_label = folder_name;
            end
            return
        end
    end
end
end

function compute_presynaptic_maps(session_paths, tt_file)

ensure_cluster_cutting_paths();

session_paths = session_paths(:);
for idx = 1:numel(session_paths)
    session_paths{idx} = ensure_trailing_filesep(session_paths{idx});
end

equalPlot('ccg_review_temp.txt', 0, ...
    'sessions', session_paths, ...
    'F', {tt_file}, ...
    'trackWidthCm', 50, ...
    'trackHeightCm', 50, ...
    'foldertag', '_ccgReview', ...
    'img_text', 'off', ...
    'suppressRedraw', 1, ...
    'exportMode', 'fast');
end

function ensure_cluster_cutting_paths()

persistent cluster_paths_ready

if isempty(cluster_paths_ready) || exist('equalPlot', 'file') ~= 2 || ...
        exist('MClustLoadVideoTrackerNvtMatlab', 'file') ~= 2
    config = load_classification_config('');
    add_dependency_path(config.mclustPath, 'MClust');
    cluster_paths_ready = true;
end

if exist('equalPlot', 'file') ~= 2
    error('equalPlot was not found. Set mclustPath in classification_config.json or pass export_options.mclustPath.');
end
end

function path_out = ensure_trailing_filesep(path_in)

path_out = char(path_in);
if isempty(path_out)
    return
end
if path_out(end) ~= filesep
    path_out = [path_out, filesep];
end
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
