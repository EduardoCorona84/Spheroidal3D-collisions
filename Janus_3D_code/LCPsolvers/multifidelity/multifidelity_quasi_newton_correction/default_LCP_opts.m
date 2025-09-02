function [opts, info] = default_LCP_opts(opts, x0, multi)
    n = numel(x0);

    if nargin < 3
        multi = false;
    end

    if ~exist('opts','var') || isempty(opts)
        opts = struct();
    end

    if ~isfield(opts, 'solver')
        opts.solver = 'bbpgd';
    end

    if ~isfield(opts, 'n')
        opts.n = n;
    end

    if ~isfield(opts, 'max_iter') && multi == false
        opts.max_iter = 100;
    end

    if ~isfield(opts, 'tol_rel')
        opts.tol_rel = 1e-4;
    end

    if ~isfield(opts, 'tol_abs')
        opts.tol_abs = 1e-6;
    end

    if ~isfield(opts, 'gamma')
        opts.gamma = 0.8;
    end

    if ~isfield(opts, 'tau')
        % warning("Tau should be set as the 1/Lhat the the best estimate of the lipschitz constant of the the A mat");
        opts.tau = 1;
    end

    if ~isfield(opts, 'tau_min')
        opts.tau_min = 1e-14;
    end

    if ~isfield(opts, 'tau_max')
        opts.tau_max = Inf;
    end

    if ~isfield(opts, 'stepSize')
        opts.stepSize = struct('init', 'uniform',...
            'eta', 'opt');
        switch lower(opts.solver)
            case 'bbpgd'
                opts.stepSize.kappa = 'bb1';
            otherwise 
                opts.stepSize.kappa = 'uniform';
        end
    end

    if ~isfield(opts, 'r')
        opts.r = min(20, n);
    else 
        assert(opts.r > 0, 'Memory/effective-rank or hessian must by positive');
        opts.r = min(opts.r, n);
    end

    opts.S = zeros(n, opts.r);
    opts.Y = zeros(n, opts.r);

    if ~isfield(opts, 'qnUpdate') || isempty(opts.qnUpdate)
        opts.qnUpdate = 'bfgs';
    else 
        assert( strcmpi('bfgs', opts.qnUpdate) ...
            || strcmpi('sr1', opts.qnUpdate), 'Must choose SR1 or BFGS')
    end

    
    if ~isfield(opts, 'prox')|| isempty(opts.prox)
        opts.prox = struct( ...
            'maxiter', 1000, ...
            'res_abstol', 0, ...
            'res_reltol', 0, ...
            'alp_abstol', 0, ...
            'alp_reltol', 1e-12, ...
            'verbose', false, ...
            'runCVX', false);
    end


%% initialize info

    if multi == true
        %do nothing, we'll deal with storage in the multi fidelity set default opts.
        info = [];
        return
    else
        info = struct('kkt', [], ...
            'iter',[],...
            'flag', [],...
            'msg',[]);
        if ~isfield(opts, 'errFcn') || isempty(opts.errFcn)
            opts.errFcn = [];
        elseif isa(opts.errFcn,'function_handle')
            info.errHist = zeros(opts.max_iter+1,1);
            info.errHist(1) = opts.errFcn(x0);
        else
            assert(iscell(opts.errFcn) && ...
                    all(cellfun(@(f) isa(f,'function_handle'), opts.errFcn)));
            info.errHist = zeros(opts.max_iter+1,numel(opts.errFcn));
            for i = 1:numel(opts.errFcn)
                fcn = opts.errFcn{i};
                info.errHist(1,i) = fcn(x0);
            end
        end

        if ~isfield(opts, 'storeIts') || isempty(opts.storeIts)
            opts.storeIts = false;
        elseif opts.storeIts
            info.iterHist = zeros(opts.max_iter+1, n);
        end
    end
end % defaultLCPOpts
