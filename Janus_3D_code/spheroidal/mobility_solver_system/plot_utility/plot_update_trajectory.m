function plot_state = plot_update_trajectory(plot_state, Ct)
    plot_state.traj_counter = plot_state.traj_counter + 1;

    % Update every trajectory for each spheroid
    for k=1:plot_state.n3
        if isempty(plot_state.traj_history{k})
            plot_state.traj_history{k} = Ct(k,:);
        else
            plot_state.traj_history{k} = [plot_state.traj_history{k} ; Ct(k,:)];
        end

        % Keep the previous max_points trajectory points
        if ~isempty(plot_state.traj_max_points) && size(plot_state.traj_history{k},1) > plot_state.traj_max_points
            plot_state.traj_history{k} = plot_state.traj_history{k}(end - plot_state.traj_max_points + 1:end, :);
        end
    end
end