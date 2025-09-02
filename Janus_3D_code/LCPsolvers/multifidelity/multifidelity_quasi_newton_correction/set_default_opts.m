function [opts, info] = set_default_opts(opts, x0)
    %this function sets default opts for the solver 
    n = numel(x0);

    if ~isfield(opts, 'outer')
        opts.outer = struct();
    end

    if ~isfield(opts.outer, 'max_iter')
        opts.outer.max_iter = 200;
    end

    if ~isfield(opts.outer, 'correction')
        opts.outer.correction = true;
    end

    %make correction opts struct
    if opts.outer.correction == true 

        if ~isfield(opts.outer, 'correction_opts')
            opts.outer.correction_opts = struct();
        end

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

        %make storage matrices for the different types of memory

        switch(opts.outer.correction_opts.memory)
            case 'dense full'
                %no r needed, just a dense matrix
                opts.outer.correction_opts.update_matrix = zeros(n, n);

            case 'dense limited'
                if ~isfield(opts.outer.correction_opts, 'r')
                    opts.outer.correction_opts.r = min(20, n);
                end

                if ~isfield(opts.outer.correction_opts, 'S')
                    opts.outer.correction_opts.S = zeros(n, opts.outer.correction_opts.r);
                end

                if ~isfield(opts.outer.correction_opts, 'Y')
                    opts.outer.correction_opts.Y = zeros(n, opts.outer.correction_opts.r);
                end

            case 'compact'
                if ~isfield(opts.outer.correction_opts, 'r')
                    opts.outer.correction_opts.r = min(20, n);
                end

                if ~isfield(opts.outer.correction_opts, 'S')
                    opts.outer.correction_opts.S = zeros(n, opts.outer.correction_opts.r);
                end

                if ~isfield(opts.outer.correction_opts, 'Y')
                    opts.outer.correction_opts.Y = zeros(n, opts.outer.correction_opts.r);
                end
                
            otherwise
                error('Unknown memory option');
        end

    end

    if ~isfield(opts.outer, 'store_updates')
        opts.outer.store_updates = false;
    end


    if ~isfield(opts.outer, 'adaptive')
        opts.outer.adaptive = 'high';
    end

    if ~isfield(opts.outer, 'solver_opts')
        %default to prox for the outer solver
        opts.outer.solver_opts = struct();
        opts.outer.solver_opts.solver = 'prox';
        %we can pass this to standard default opts with a flat that will tell it to exclude certain things
        [info.outer, opts.outer.solver_opts] = default_LCP_opts(opts.outer.solver_opts, x0, true);
    end

    %now all the outer solver opts should be set, we want to set the info struct for outer iterates ourselves

    if ~isfield(opts.outer, 'storeIts') || isempty(opts.outer.storeIts)
        opts.outer.storeIts = false;
    elseif opts.outer.storeIts
        info.outer.iterHist = zeros(opts.outer.max_iter+1, n);
    end

    %add in other outer info struct later

    %now we go through inner opts
    if ~isfield(opts.inner, 'enabled')
        opts.inner.enabled = true;
    end


    if ~isfield(opts.inner, 'solver')
        opts.inner.solver = 'outer preconditioned prox';
        %we are going to use the exact same information for the prox steps as the outer solver. some of the other things might have to be edited a bit.
        opts.inner.solver_opts = opts.outer.solver_opts;
    end
    
    if ~isfield(opts.inner, 'max_iter')
        opts.inner.solver_opts.max_iter = 1;
    end

    %for every outer iteration, we have an entry in the cell array to store the info the inner iterations corresponding to the outer iteration.
    info.inner = cell(opts.outer.max_iter + 1, 1);

end