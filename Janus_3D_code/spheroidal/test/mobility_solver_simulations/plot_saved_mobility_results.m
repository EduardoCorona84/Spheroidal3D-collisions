function replay_info = plot_saved_mobility_results(source_file_loc, opts)
%{
Replay and summarize saved spheroidal mobility runs.

Usage:
    plot_saved_mobility_results('path/to/results.mat')
    plot_saved_mobility_results('path/to/results.mat', struct('plotColor','sigma'))

Inputs:
    source_file_loc
        - file path to a .mat file containing tt, Xt, and Ct.
    opts (optional)
        Replay controls:
            replay          (default true)
            frameStep       (default 1)
            replayStartStep (default 0)
            replayEndStep   (default inf)
            pauseTime       (default 0.02)
            finalPauseTime (default 0; also applied to saved animation)
            saveAnimation   (default false)
            animationFormat ('gif',  'mp4', 'gifv', default 'gif')
            animationPath   (default auto name in pwd)
            frameRate       (default auto from pauseTime)
            gifLoopCount    (default inf)
            gifDelayTime    (default 1/frameRate)
        Plot controls:
            See plot_init.

Output:
    replay_info struct with selected file and frame indices.
%}

if nargin < 1 || isempty(source_file_loc)
    error('source_file_loc is required. Pass a file path to a .mat mobility checkpoint.');
end
if nargin < 2 || isempty(opts)
    opts = struct();
end

LOCAL_add_required_paths();

settings = LOCAL_parse_settings(opts);
plot_params = LOCAL_parse_plot_params(opts);

result_file = LOCAL_validate_result_file(source_file_loc);
snapshot = LOCAL_load_data(result_file);

state_indices = LOCAL_filter_state_indices(snapshot.state_indices, settings);
data = snapshot.data;
n3 = snapshot.n3;
np = snapshot.np;

if isempty(state_indices)
    error('No timesteps fall within replayStartStep/replayEndStep for the selected file.');
end

if settings.replay
    plot_state = plot_init(plot_params, data.Xt{state_indices(1)}, np, n3);
    if plot_state.traj_enable
        plot_state = plot_update_trajectory(plot_state, data.Ct{state_indices(1)});
    end

    writer = LOCAL_init_animation_writer(settings);
    frame_indices = 1:settings.frameStep:numel(state_indices);
    for frame_id = frame_indices
        state_idx = state_indices(frame_id);
        t = LOCAL_get_time(data.tt, state_idx);
        [sigma_k, mu_k, FT_k] = LOCAL_get_plot_data(data, state_idx);

        if plot_state.traj_enable && frame_id > 1
            plot_state = plot_update_trajectory(plot_state, data.Ct{state_idx});
        end
        plot_state = plot_update(plot_state, data.Xt{state_idx}, data.Ct{state_idx}, t, sigma_k, mu_k, FT_k);
        if plot_params.plotShowTimestep
            LOCAL_update_title_with_timestep(plot_state, state_idx, t);
        end
        writer = LOCAL_write_animation_frame(writer, plot_state.fig);

        if settings.pauseTime > 0 && frame_id < frame_indices(end)
            pause(settings.pauseTime);
        end
    end
    if settings.finalPauseTime > 0
        pause(settings.finalPauseTime);
    end
    writer = LOCAL_close_animation_writer(writer);
else
    writer = LOCAL_init_animation_writer(settings);
    if writer.enabled
        warning('saveAnimation=true ignored because replay=false.');
    end
end

replay_info = struct( ...
    'file', snapshot.file, ...
    'state_indices', state_indices, ...
    'n3', n3, ...
    'np', np, ...
    'final_time', snapshot.last_time, ...
    'animation_file', writer.output_path ...
);

fprintf('\nLoaded %s', snapshot.file);
fprintf('\nFrames replayed: %d', numel(state_indices));
fprintf('\nFinal time: %.6g\n', snapshot.last_time);
if ~isempty(writer.output_path)
    fprintf('Animation saved: %s\n', writer.output_path);
end
end %% END MAIN FUNCTION

function settings = LOCAL_parse_settings(opts)
    % Setup default struct
    settings = struct( ...
        'replay', true, ...
        'frameStep', 1, ...
        'replayStartStep', 0, ...
        'replayEndStep', inf, ...
        'pauseTime', 0.02, ...
        'finalPauseTime', 0, ...
        'saveAnimation', false, ...
        'animationFormat', 'gif', ...
        'animationPath', '', ...
        'frameRate', [], ...
        'gifLoopCount', inf, ...
        'gifDelayTime', [] ...
    );

    names = fieldnames(settings);
    for i = 1:numel(names)
        key = names{i};
        if isfield(opts, key)
            settings.(key) = opts.(key);
        end
    end

    % Normalize/validate text settings.
    if isstring(settings.animationFormat)
        if ~isscalar(settings.animationFormat)
            error('animationFormat must be a single text value.');
        end
        settings.animationFormat = char(settings.animationFormat);
    end
    if ~ischar(settings.animationFormat)
        error('animationFormat must be text.');
    end
    settings.animationFormat = lower(strtrim(settings.animationFormat));
    if isempty(settings.animationFormat)
        settings.animationFormat = 'gif';
    end

    if isstring(settings.animationPath)
        settings.animationPath = char(settings.animationPath);
    end
    if ~ischar(settings.animationPath)
        error('animationPath must be text.');
    end

    if isempty(settings.frameRate)
        if isscalar(settings.pauseTime) && settings.pauseTime > 0
            settings.frameRate = 1 / settings.pauseTime;
        else
            settings.frameRate = 30;
        end
    end

    if isempty(settings.gifDelayTime)
        settings.gifDelayTime = 1 / settings.frameRate;
    end
end

function plot_params = LOCAL_parse_plot_params(opts)
    % Default struct
    plot_params = struct( ...
        'plotFlag', true, ...
        'plotTrajectories', true, ...
        'plotTrajMaxPoints', 200, ...
        'plotSurfaceAlpha', 0.2, ...
        'plotGrid', true, ...
        'plotColor', 'mu', ...
        'plotColorMode', 'inf', ...
        'plotColormap', 'parula', ...
        'plotForceVectors', true, ...
        'plotView', [45 0], ...
        'plotShowTimestep', false ...
    );

    fn = fieldnames(opts);
    for i = 1:numel(fn)
        key = fn{i};
        if startsWith(key, 'plot')
            plot_params.(key) = opts.(key);
        end
    end
end

function state_indices = LOCAL_filter_state_indices(all_state_indices, settings)
    % Extract the requested timesteps
    state_steps = max(0, all_state_indices - 1);
    keep = state_steps >= settings.replayStartStep;
    if isfinite(settings.replayEndStep)
        keep = keep & (state_steps <= settings.replayEndStep);
    end
    state_indices = all_state_indices(keep);
end

function LOCAL_update_title_with_timestep(plot_state, state_idx, t)
    step_idx = max(0, state_idx - 1);

    % Prevent crashing when closing plot window
    if ~isfield(plot_state, 'ax') || ~ishandle(plot_state.ax)
        return;
    end

    if isfield(plot_state, 'title_prefix') && ~isempty(plot_state.title_prefix)
        title(plot_state.ax, sprintf('%s step %d, t = %.4f', plot_state.title_prefix, step_idx, t));
    else
        title(plot_state.ax, sprintf('step %d, t = %.4f', step_idx, t));
    end
end

function LOCAL_add_required_paths()
    this_file = mfilename('fullpath');
    this_dir = fileparts(this_file);
    repo_root = fileparts(fileparts(fileparts(this_dir)));

    addpath(genpath(fullfile(repo_root, 'spheroidal')));
    addpath(genpath(fullfile(repo_root, 'support')));
end

function result_loc = LOCAL_validate_result_file(source_file_loc)
    source_file_loc = char(source_file_loc);

    result_loc = strtrim(source_file_loc);
    if ~exist(result_loc, 'file')
        error('source_file_loc file does not exist: %s', result_loc);
    end

    if isempty(result_loc)
        error('source_file_loc cannot be empty.');
    end

    [~, ~, ext] = fileparts(result_loc);
    if ~strcmpi(ext, '.mat')
        error('source_file_loc must point to a .mat file: %s', result_loc);
    end
end

function data = LOCAL_load_data(result_file)
    S = load(result_file);
    ok = isfield(S, 'tt') && isfield(S, 'Xt') && isfield(S, 'Ct');
    if ~ok
        error('Selected file is missing required fields tt, Xt, Ct: %s', result_file);
    end

    [state_indices, n3, np] = LOCAL_valid_state_indices(S);
    if isempty(state_indices)
        error('Selected file has no valid mobility states to replay: %s', result_file);
    end

    last_idx = state_indices(end);
    last_time = LOCAL_get_time(S.tt, last_idx);
    data = struct( ...
        'file', result_file, ...
        'data', S, ...
        'state_indices', state_indices, ...
        'last_time', last_time, ...
        'last_idx', last_idx, ...
        'n3', n3, ...
        'np', np ...
    );
end

function [state_indices, n3, np] = LOCAL_valid_state_indices(S)
    nt = numel(S.tt);
    nx = numel(S.Xt);
    nc = numel(S.Ct);
    n = min([nt nx nc]);

    valid = false(n, 1);
    n3 = 0;
    np = 0;

    for i = 1:n
        Xt_i = S.Xt{i};
        Ct_i = S.Ct{i};

        if isempty(Xt_i) || isempty(Ct_i)
            continue;
        end
        if size(Ct_i, 2) ~= 3 || size(Xt_i, 2) ~= 3
            continue;
        end

        n3_i = size(Ct_i, 1);
        if n3_i <= 0
            continue;
        end
        np_i = size(Xt_i, 1) / n3_i;
        if abs(np_i - round(np_i)) > 1e-12
            continue;
        end

        valid(i) = true;
        n3 = n3_i;
        np = round(np_i);
    end

    state_indices = find(valid);
end

function t = LOCAL_get_time(tt, idx)
    if isempty(tt)
        t = idx - 1;
        return;
    end

    tt = tt(:);
    if idx <= numel(tt) && isfinite(tt(idx))
        t = tt(idx);
    else
        t = idx - 1;
    end
end

function [sigma_k, mu_k, FT_k] = LOCAL_get_plot_data(S, state_idx)
    if state_idx <= 1
        sigma_k = [];
        mu_k = [];
        FT_k = [];
        return;
    end

    step_idx = state_idx - 1;
    sigma_k = LOCAL_get_cell(S, 'sigma', step_idx);
    mu_k = LOCAL_get_cell(S, 'mu', step_idx);
    FT_k = LOCAL_get_cell(S, 'FT', step_idx);
end

function val = LOCAL_get_cell(S, field_name, idx)
    val = S.(field_name){idx};
end

function writer = LOCAL_init_animation_writer(settings)
    writer = struct( ...
        'enabled', settings.saveAnimation, ...
        'mode', '', ...
        'frame_id', 0, ...
        'video', [], ...
        'last_frame', [], ...
        'path_mp4', '', ...
        'path_gif', '', ...
        'path_gifv', '', ...
        'output_path', '', ...
        'gif_loop_count', settings.gifLoopCount, ...
        'gif_delay_time', settings.gifDelayTime, ...
        'final_pause_time', settings.finalPauseTime, ...
        'last_imind', [], ...
        'last_cmap', [] ...
    );

    if ~writer.enabled
        return;
    end

    [path_mp4, path_gif, path_gifv, mode, output_path] = ...
        LOCAL_resolve_animation_paths(settings.animationPath, settings.animationFormat);
    writer.mode = mode;
    writer.path_mp4 = path_mp4;
    writer.path_gif = path_gif;
    writer.path_gifv = path_gifv;
    writer.output_path = output_path;

    if strcmp(mode, 'video')
        writer.video = VideoWriter(writer.path_mp4, 'MPEG-4');
        writer.video.FrameRate = settings.frameRate;
        open(writer.video);
    end
end

function writer = LOCAL_write_animation_frame(writer, fig)
    if ~writer.enabled || ~ishandle(fig) % Prevent crash
        return;
    end

    frame = getframe(fig);
    writer.frame_id = writer.frame_id + 1;
    writer.last_frame = frame;

    switch writer.mode
        case 'video'
            writeVideo(writer.video, frame);
        case 'gif' % Ripped from 'imwrite' docs
            img = frame2im(frame);
            [imind, cmap] = rgb2ind(img, 256);
            writer.last_imind = imind;
            writer.last_cmap = cmap;
            if writer.frame_id == 1
                imwrite(imind, cmap, writer.path_gif, 'gif', ...
                    'LoopCount', writer.gif_loop_count, ...
                    'DelayTime', writer.gif_delay_time);
            else
                imwrite(imind, cmap, writer.path_gif, 'gif', ...
                    'WriteMode', 'append', ...
                    'DelayTime', writer.gif_delay_time);
            end
        otherwise
            error('Unknown animation mode "%s".', writer.mode);
    end
end

function writer = LOCAL_close_animation_writer(writer)
    if writer.final_pause_time > 0 && writer.frame_id > 0
        writer = LOCAL_append_final_frame_pause(writer);
    end

    if strcmp(writer.mode, 'video') && ~isempty(writer.video)
        close(writer.video);
    end

    if strcmp(writer.mode, 'video') && ~isempty(writer.path_gifv)
        copyfile(writer.path_mp4, writer.path_gifv, 'f');
    end
end

function writer = LOCAL_append_final_frame_pause(writer)
    % Add final frame to the movie for a pause effect
    switch writer.mode
        case 'video'
            if isempty(writer.video) || isempty(writer.last_frame)
                return;
            end
            extra_frames = max(0, round(writer.final_pause_time * writer.video.FrameRate));
            for k = 1:extra_frames
                writeVideo(writer.video, writer.last_frame);
            end
        case 'gif'
            if isempty(writer.last_imind) || isempty(writer.last_cmap)
                return;
            end
            base_delay = max(eps, writer.gif_delay_time);
            extra = writer.final_pause_time;
            n_full = floor(extra / base_delay);
            rem_delay = extra - n_full*base_delay;
            for k = 1:n_full
                imwrite(writer.last_imind, writer.last_cmap, writer.path_gif, 'gif', ...
                    'WriteMode', 'append', ...
                    'DelayTime', base_delay);
            end
            if rem_delay > 1e-12
                imwrite(writer.last_imind, writer.last_cmap, writer.path_gif, 'gif', ...
                    'WriteMode', 'append', ...
                    'DelayTime', rem_delay);
            end
        otherwise
            error('Invalid mode given.');
    end
end

function [path_mp4, path_gif, path_gifv, mode, output_path] = LOCAL_resolve_animation_paths(animation_path, animation_format)
    path_mp4 = '';
    path_gif = '';
    path_gifv = '';

    if isempty(animation_format)
        fmt = 'gif';
    else
        fmt = strtrim(animation_format);
    end

    if ~isempty(animation_path)
        [folder, base_name, ext] = fileparts(animation_path);
        if isempty(folder), folder = '.'; end
        base_path = fullfile(folder, base_name);
    else
        error('Animation path not given.');
    end

    switch fmt
        case 'gif'
            mode = 'gif';
            if isempty(ext)
                path_gif = [base_path '.gif'];
            else
                path_gif = [base_path ext];
            end
            output_path = path_gif;
        case 'mp4'
            mode = 'video';
            if isempty(ext)
                path_mp4 = [base_path '.mp4'];
            else
                path_mp4 = [base_path ext];
            end
            output_path = path_mp4;
        case 'gifv'
            mode = 'video';
            path_mp4 = [base_path '.mp4'];
            path_gifv = [base_path '.gifv'];
            output_path = path_gifv;
        otherwise
            error('Unsupported animationFormat "%s". Use gif, mp4, or gifv.', fmt);
    end
end
