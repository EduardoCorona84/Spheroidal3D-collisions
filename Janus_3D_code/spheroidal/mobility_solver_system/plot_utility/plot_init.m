function plot_state = plot_init(Fparams, Xt, np, n3)
    plot_state = struct();
    plot_state.fig = figure;
    plot_state.ax = axes('Parent', plot_state.fig);
    plot_state.np = np;
    plot_state.n3 = n3;

    plot_state.view = [-1 0.5 0.5];
    if isfield(Fparams,'plotView')
        plot_state.view = Fparams.plotView.';
    end

    plot_state.view_init = [];
    plot_state.view_applied = false;
    if isfield(Fparams,'plotViewInit')
        plot_state.view_init = Fparams.plotViewInit.';
    end

    plot_state.axis = [];
    if isfield(Fparams,'plotAxis')
        plot_state.axis = Fparams.plotAxis.';
    end

    mins = min(Xt, [], 1);
    maxs = max(Xt, [], 1);
    plot_state.diam = max(maxs - mins);
    if plot_state.diam <= 0
        plot_state.diam = 1;
    end

    plot_state.title_prefix = '';
    if isfield(Fparams,'plotTitle')
        plot_state.title_prefix = Fparams.plotTitle;
    end

    plot_state.color_source = 'none';
    if isfield(Fparams,'plotColor')
        plot_state.color_source = Fparams.plotColor;
    end
    
    plot_state.color_mode = 'inf';
    if isfield(Fparams,'plotColorMode')
        plot_state.color_mode = Fparams.plotColorMode;
    end

    plot_state.color_limits = [];
    if isfield(Fparams,'plotColorLimits')
        plot_state.color_limits = Fparams.plotColorLimits.';
    end

    plot_state.traj_enable = false;
    if isfield(Fparams,'plotTrajectories')
        plot_state.traj_enable = Fparams.plotTrajectories;
    end

    plot_state.traj_max_points = [];
    if isfield(Fparams,'plotTrajMaxPoints')
        plot_state.traj_max_points = max(200, Fparams.plotTrajMaxPoints);
    end

    plot_state.traj_line_width = 1.0;

    plot_state.traj_counter = 0;
    plot_state.traj_history = cell(n3,1);
    plot_state.traj_lines = gobjects(n3,1);
    plot_state.traj_colors = lines(n3);

    plot_state.surface_alpha = [];
    if isfield(Fparams,'plotSurfaceAlpha')
        plot_state.surface_alpha = Fparams.plotSurfaceAlpha;
    end

    plot_state.grid_on = true;
    if isfield(Fparams,'plotGrid')
        plot_state.grid_on = Fparams.plotGrid;
    end
end