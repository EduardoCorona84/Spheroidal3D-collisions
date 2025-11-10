function [opts, info] = set_default_opts(opts, x0, fg, fg_low)
    %this function sets default opts for the solver 

    %in the case in which A or b are not given, we can get them from the function handle.
    %This requires adding the path to this function
    addpath('../../utilities/');
    n = numel(x0);

    %warm start opts
    if ~isfield(opts, 'warm')
        opts.warm = struct();
    end

    if ~isfield(opts.warm, 'enabled')
        opts.warm.enabled = false;
    end

    if ~isfield(opts.warm, 'solver_opts')
        opts.warm.solver_opts = struct();
    end

    if ~isfield(opts.warm.solver_opts, 'type')
        opts.warm.solver_opts.type = 'bbpgd';
    end

    if ~isfield(opts.warm.solver_opts, 'max_iter')
        opts.warm.solver_opts.max_iter = 20;
    end
    
    %we won't make the info struct as the warm start solver will do this automatically. We will just assign the resulting struct to the info struct of the multifielity (under warm).
    if ~isfield(opts.warm.solver_opts, 'storeIts') || isempty(opts.warm.solver_opts.storeIts)
        opts.warm.solver_opts.storeIts = true;
    end


    %outer solver opts
    if ~isfield(opts, 'outer')
        opts.outer = struct();
    end

    if ~isfield(opts.outer, 'max_iter')
        opts.outer.max_iter = 50;
    end

    if ~isfield(opts.outer, 'correction')
        opts.outer.correction = true;
    end

    %if we warm start, we can use this as a correction on the first iteration of multifidelity solve
    if opts.warm.enabled == true && ~isfield(opts.outer, 'warm_correction')
        opts.outer.warm_correction = true;
    end

    %warm correction is used for if correcting off a warm start or correcting off a nonzero starting point
    if opts.warm.enabled == false && ~isfield(opts.outer, 'warm_correction')
        opts.outer.warm_correction = false;
    end

    %here we need to check if we have x0 = 0 for the multifidelity part. 
    %If x0 is zero and opts.outer.warm_correction == true, the update will fail
    if all(x0 == 0) && opts.outer.warm_correction == true
        opts.outer.warm_correction = false;
    end

    %now add logic for correcting along b
    if opts.warm.enabled == true && ~isfield(opts.outer, 'b_correction')
        opts.outer.b_correction = true;
    end

    if opts.warm.enabled == false && ~isfield(opts.outer, 'b_correction')
        opts.outer.b_correction = false;
    end

    %in the case with no warm start and x0 equals 0, b correction is redundant.
    if (opts.warm.enabled == false && all(x0 == 0)) && opts.outer.b_correction == true 
        opts.outer.b_correction = false;
    end


    %these are out here as right now I am using the same secant directions for correction and prox. This could be changed in the future.

    if ~isfield(opts.outer, 'correction')
        opts.outer.correction = true;
    end
    
    if ~isfield(opts.outer, 'correction_opts')
        opts.outer.correction_opts = struct();
    end

    if ~isfield(opts.outer.correction_opts, 'direction')
        opts.outer.correction_opts.direction = 'secant';
    end



    %make fill in the correction opts struct
    if opts.outer.correction == true 


        if ~isfield(opts.outer.correction_opts, 'update')
            opts.outer.correction_opts.update = 'sr1';
        end
         
        %no skips by default
        if ~isfield(opts.outer.correction_opts, 'skip')
            opts.outer.correction_opts.skips = false;
        end

        if ~isfield(opts.outer.correction_opts, 'memory')
            opts.outer.correction_opts.memory = 'dense full';
        end

        %make storage matrices for the different types of memory, also depends on the update type

        switch(opts.outer.correction_opts.memory)
            case 'dense full'
                %no r needed, just a dense matrix
                opts.outer.correction_opts.update_matrix = zeros(n, n);

            case 'dense limited'
                if ~isfield(opts.outer.correction_opts, 'r')
                    opts.outer.correction_opts.r = min(20, n);
                end

                opts.outer.correction_opts.curr_mem = 0;

                if ~isfield(opts.outer.correction_opts, 'S')
                    opts.outer.correction_opts.S = zeros(n, opts.outer.correction_opts.r);
                end

                if ~isfield(opts.outer.correction_opts, 'Y')
                    opts.outer.correction_opts.Y = zeros(n, opts.outer.correction_opts.r);
                end

                %add an additional flag if we want to store the first update persistently. This is in the case that we want to store along the direction b consistently.
                if ~isfield(opts.outer.correction_opts, 'persistent_first_update')
                    opts.outer.correction_opts.persistent_first_update = false;
                end

            case 'compact'
                if ~isfield(opts.outer.correction_opts, 'r')
                    opts.outer.correction_opts.r = min(20, n);
                end

                opts.outer.correction_opts.curr_mem = 0;

                if ~isfield(opts.outer.correction_opts, 'S')
                    opts.outer.correction_opts.S = zeros(n, opts.outer.correction_opts.r);
                end

                if ~isfield(opts.outer.correction_opts, 'utility_matrix')
                    opts.outer.correction_opts.utility_matrix = zeros(n, opts.outer.correction_opts.r);
                end

                %if we are using BFGS or DFP, we need to store an additional matrix for the compact representation (kinda, should go over this)
                if strcmpi(opts.outer.correction_opts.update, 'dfp') || strcmpi(opts.outer.correction_opts.update, 'bfgs')
                    opts.outer.correction_opts.Y = zeros(n, opts.outer.correction_opts.r);
                end

                %add an additional flag if we want to store the first update persistently
                if ~isfield(opts.outer.correction_opts, 'persistent_first_update')
                    opts.outer.correction_opts.persistent_first_update = false;
                end

            otherwise
                error('Unknown memory option');
        end

    end

    if ~isfield(opts.outer, 'store_updates')
        opts.outer.store_updates = false;
    end


    if ~isfield(opts.outer, 'adaptive')
        opts.outer.adaptive = 'retry';
    end

    if ~isfield(opts.outer, 'solver_opts')
        %default to prox for the outer solver
        opts.outer.solver_opts = struct();
    end

    if ~isfield(opts.outer.solver_opts, 'solver')
        opts.outer.solver_opts.solver = 'prox';
    end

    %we can pass this to standard default opts with a flat that will tell it to exclude certain things
    [opts.outer.solver_opts, ~] = default_LCP_opts(opts.outer.solver_opts, x0, true);

    %if taking optimal step sizes, add the necessary info to the opts struct
    if strcmpi(opts.outer.solver_opts.stepSize.eta, 'opt')
        % If we are using the optimal step size, we need to set the
        % appropriate options in the solver_opts
        if ~isfield(opts.outer.solver_opts, 'A')
            opts.outer.solver_opts.A = create_A(fg);
        end

        if ~isfield(opts.outer.solver_opts, 'b')
            %get b from a function evaluation, this is really not something that should be done as a function eval can be very expensive, this is used as a fallback.
            [~, grad_ones, A_ones] = fg(ones(n,1), [], [], []);
            opts.outer.solver_opts.b = grad_ones - A_ones;
        end

    end

    %now all the outer solver opts should be set, we want to set the info struct for outer iterates ourselves 

    if ~isfield(opts.outer, 'storeIts') || isempty(opts.outer.storeIts)
        opts.outer.storeIts = true;
    end

    if ~isfield(opts.outer, 'storeFuncs')
        opts.outer.storeFuncs = false;
    end

    if opts.outer.storeIts == true
        info.outer.iterHist = zeros(opts.outer.max_iter+1, n);
    end

    if opts.outer.storeFuncs == true
        info.outer.funcHist = zeros(opts.outer.max_iter+1, 1);
    end

    if opts.outer.storeKKT == true
        info.outer.kktHist = zeros(opts.outer.max_iter+1, 1);
    end

    %add in other outer info struct later


    %now we go through inner opts
    if ~isfield(opts.inner, 'enabled')
        opts.inner.enabled = true;
    end

    %by default us outer preconditioned prox, we won't set all the opts as the function does it for us (in outer_solver)
    if ~isfield(opts.inner, 'solver_opts')
        opts.inner.solver_opts = struct();
    end

    if ~isfield(opts.inner.solver_opts, 'solver')
        opts.inner.solver_opts.solver = 'outer preconditioned prox';
    end

    if ~isfield(opts.inner.solver_opts, 'A')
        opts.inner.solver_opts.A = create_A(fg_low);
    end

    if ~isfield(opts.inner.solver_opts, 'max_iter')
        opts.inner.solver_opts.max_iter = 1;
    end

    if ~isfield(opts.inner.solver_opts, 'gradient_mode')
        opts.inner.solver_opts.gradient_mode = 'full'; %can also be 'step' or 'reuse'
    end

    if ~isfield(opts.inner.solver_opts, 'adaptive')
        opts.inner.solver_opts.adaptive = struct();
    end

    if ~isfield(opts.inner.solver_opts.adaptive, 'enabled')
        opts.inner.solver_opts.adaptive.enabled = false;
    end

    if ~isfield(opts.inner.solver_opts.adaptive, 'abs_err_estimate')
        opts.inner.solver_opts.adaptive.abs_err_estimate = 0.2;
    end

    if ~isfield(opts.inner.solver_opts.adaptive, 'descent_parameter')
        opts.inner.solver_opts.adaptive.descent_parameter = 1/100;
    end

    opts.inner.solver_opts.b = opts.outer.solver_opts.b;

    opts.warm.solver_opts.A = opts.inner.solver_opts.A;

    opts.warm.solver_opts.b = opts.inner.solver_opts.b;

    %for every outer iteration, we have an entry in the cell array to store the info the inner iterations corresponding to the outer iteration.
    info.inner = cell(opts.outer.max_iter + 1, 1);

end