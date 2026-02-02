function plot_state = plot_update(plot_state, Xt, Ct, t, sigma, mu)
    % Below is needed to stop crashing after closing the plot window.
    if ~ishandle(plot_state.fig) || ~ishandle(plot_state.ax)
        return;
    end

    Xin = LOCAL_plot_points(Xt, Ct, plot_state.np, plot_state.n3);
    axes(plot_state.ax);
    hold(plot_state.ax, 'off');
    C = [];

    switch plot_state.color_source
        case 'sigma'
            C = LOCAL_plot_colors(sigma, plot_state.np, plot_state.n3, plot_state.color_mode);
        case 'mu'
            C = LOCAL_plot_colors(mu, plot_state.np, plot_state.n3, plot_state.color_mode);
    end

    if isempty(C)
        h = plotb(Xin);
    else
        h = plotb(Xin, C);
        if ~isempty(plot_state.color_limits)
            clim(plot_state.ax, plot_state.color_limits);
        end
    end

    if ~isempty(plot_state.surface_alpha)
        set(h, 'FaceAlpha', plot_state.surface_alpha);
    end

    if plot_state.traj_enable
        plot_state = LOCAL_plot_trajectories(plot_state, Ct);
    end

    if isempty(plot_state.axis)
        axis(plot_state.ax, LOCAL_plot_axis(Ct, plot_state.diam));
    else
        axis(plot_state.ax, plot_state.axis);
    end

    if ~plot_state.view_applied && ~isempty(plot_state.view_init)
        view(plot_state.ax, plot_state.view_init);
        plot_state.view_applied = true;
    elseif isempty(plot_state.view_init) && ~isempty(plot_state.view)
        view(plot_state.ax, plot_state.view);
    end

    grid(plot_state.ax, plot_state.grid_on);

    if ~isempty(plot_state.title_prefix)
        title(plot_state.ax, sprintf('%s t = %.4f', plot_state.title_prefix, t));
    else
        title(plot_state.ax, sprintf('t = %.4f', t));
    end
    drawnow;
end

function Xin = LOCAL_plot_points(Xt, Ct, np, n3)
    Xin = zeros(3*np, n3);
    for k=1:n3
        idx = (1:np) + np*(k-1);
        Xk = Xt(idx,:) + Ct(k,:);
        Xin(:,k) = [Xk(:,1); Xk(:,2); Xk(:,3)];
    end
end

function C = LOCAL_plot_colors(data, np, n3, mode)
    C = [];
    if isempty(data)
        return;
    end

    data3 = reshape(data, 3, []);
    switch mode
        case 'l2'
            mag = sqrt(sum(abs(data3).^2, 1));
        otherwise
            mag = max(abs(data3), [], 1);
    end
    C = reshape(mag, np, n3);
end

function plot_state = LOCAL_plot_trajectories(plot_state, Ct)
    if isempty(plot_state.traj_history)
        return;
    end

    for k=1:plot_state.n3
        traj = plot_state.traj_history{k};
        if size(traj,1) < 2
            continue;
        end

        hold(plot_state.ax, 'on');
        if ishandle(plot_state.traj_lines(k))
            set(plot_state.traj_lines(k), 'XData', traj(:,1), 'YData', traj(:,2), 'ZData', traj(:,3));
        else
            plot_state.traj_lines(k) = plot3(plot_state.ax, traj(:,1), traj(:,2), traj(:,3), ...
                'Color', plot_state.traj_colors(k,:), ...
                'LineWidth', plot_state.traj_line_width ...
            );
        end
    end

    hold(plot_state.ax, 'off');
end

function axvec = LOCAL_plot_axis(Ct, diam)
    if size(Ct,1) > 1
        mn = min(Ct, [], 1);
        mx = max(Ct, [], 1);
        axvec = [mn - 1.5*diam; mx + 1.5*diam];
        axvec = axvec(:).';
    else
        axvec = [Ct(1,1)-1.5*diam, Ct(1,1)+1.5*diam, ...
                 Ct(1,2)-1.5*diam, Ct(1,2)+1.5*diam, ...
                 Ct(1,3)-1.5*diam, Ct(1,3)+1.5*diam];
    end
end
