function [opts, info] = set_default_opts(x0, opts)
    %Based off what nic did for prox, the idea is this function will set the default options for the multi-fidelity solver and create the info struct.

    %set default values for outer solver/iteration

    %max high fidelity evaluations
    if ~isfield(opts, 'outer')
        opts.outer = struct();
    end

    if ~isfield(opts.outer, 'max_iter')
        opts.outer.max_iter = 100;
    end

    %warm start the high fidelity solver with low fidelity evaluations
    if ~isfield(opts.outer, 'warm_start')
        opts.outer.warm_start.enabled = false;
    end

    %first step of the high fidelity solver
    if ~isfield(opts.outer, 'step_init')
        opts.outer.step_init = 0.1;
    end

    opts.outer.update_matrix = zeros(length(x0)); 

    %use DFP as default as this directly targets the matrix (not its inverse)
    if ~isfield(opts.outer, 'low_update')
        opts.outer.low_update = 'sr1';
    end

    %default to prox qn
    if ~isfield(opts.outer, 'solver')
        opts.outer.solver = 'prox';
    end

    if strcmp(opts.outer.solver, 'prox')
        if ~isfield(opts.outer, 'prox')
            opts.outer.prox = struct();
        end

        if ~isfield(opts.outer.prox, 'r')
            opts.outer.prox.r = 10; 
        end

        %initialize needed matrices in prox
        opts.outer.prox.H = eye(length(x0));
        opts.outer.prox.S = zeros(length(x0), opts.outer.prox.r);
        opts.outer.prox.Y = zeros(length(x0), opts.outer.prox.r);


        if ~isfield(opts.outer.prox, 'gamma')
            opts.outer.prox.gamma = 0.8;
        end

        if ~isfield(opts.outer.prox, 'tau')
            opts.outer.prox.tau = 1;
        end

        if ~isfield(opts.outer.prox, 'tau_min')
            opts.outer.prox.tau_min = 1e-14;
        end

        if ~isfield(opts.outer.prox, 'tau_max')
            opts.outer.prox.tau_max = Inf;
        end

        if ~isfield(opts.outer.prox, 'qnUpdate')
            opts.outer.prox.qnUpdate = 'BFGS';
        end
    end

    if ~isfield(opts, 'inner')
        opts.inner = struct();
    end

    %enable low fidelity steps, without this its just high fidelity standard prox quasi-Newton
    if ~isfield(opts.inner, 'enabled')
        opts.inner.enabled = true;
    end

    %set default solver for the inner solver
    if ~isfield(opts.inner, 'solver')
        opts.inner.solver = 'outer prox qn';
        opts.inner.prox = opts.outer.prox;
    end

    %set default max for inner solver. This is really max number of steps with corrected gradient.
    if ~isfield(opts.inner, 'max_iter')
        opts.inner.max_iter = 5;
    end

    %history for outer solver
    if ~isfield(opts.outer, 'store_evals')
        opts.outer.store_evals = false;
    end

    if ~isfield(opts.outer, 'store_grad')
        opts.outer.store_grad = false;
    end

    if ~isfield(opts.outer, 'store_f')
        opts.outer.store_f = false;
    end

    % set up info struct/history 
    info.outer.x_iters = zeros(length(x0), opts.outer.max_iter);
    info.outer.matvecs = zeros(1, opts.outer.max_iter);
    if opts.outer.store_evals
        info.outer.evals = zeros(length(x0), opts.outer.max_iter);
    end
    if opts.outer.store_grad
        info.outer.grad_iters = zeros(length(x0), opts.outer.max_iter);
    end
    if opts.outer.store_f
        info.outer.f_iters = zeros(1, opts.outer.max_iter);
    end

    %a cell array, where each array corresponds to one outer/high fidelity evaluation. A struct of the information associated with the inner solver/low fidelity will be in each cell, some may be empty.
    info.inner = cell(1, opts.outer.max_iter);

end