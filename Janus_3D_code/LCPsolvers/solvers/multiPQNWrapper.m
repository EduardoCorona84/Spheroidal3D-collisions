function [x, info, opts] = multiPQNWrapper(A_list, b_list, x0, base_opts)

    if nargin < 4 || isempty(base_opts)
        base_opts = struct();
    end
    
    numFidelities = length(A_list);
    
    % Build from the bottom up (lowest fidelity first)
    % Start with the lowest fidelity level (no .low field = base case)
    base_opts.A = A_list{numFidelities};
    base_opts.b = b_list{numFidelities};
    base_opts.solver = 'proxquasinewton';
    base_opts = defaultLCPOpts(base_opts, x0);
    opts = base_opts;
    
    % Wrap each level going up the hierarchy
    for i = numFidelities-1:-1:1
        % Set opts 
        outer_opts.A = A_list{i};
        outer_opts.b = b_list{i};
        outer_opts.solver = 'multifidelityPQN';
        outer_opts = defaultLCPOpts(outer_opts, x0);
        outer_opts.low = opts;  % Previous level becomes the .low field
        opts = outer_opts;
    end
    
    % Define the loss function for the highest fidelity
    fg = @(x, Ax) quadraticLoss(x, opts.A, opts.b, Ax);
    
    % Call the solver
    [x, info, opts] = multifidelityProxQuasiNewton(fg, x0, opts);
end