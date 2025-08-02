function [converged, info] = checkConvergence(k, f_k, x_k, ...
    grad_k, kappa, info, opts)
persistent old_kkt
if k == 0
    old_kkt = [];
end
% Check for convergence
converged = false;
if k >= opts.max_iter
    converged = true;
elseif kappa < 10*eps
    % if step direction gets too small give up
    info.flag =  5;
    converged = true;
else
    % kkt conditions / LCP being satisified is equivalent to
    % the grad_k(i) = 0 \perp x_k(i) = 0
    phi = min(grad_k,x_k);
    kkt = 0.5*dot(phi, phi);
    if ~isempty(old_kkt) && (abs(kkt - old_kkt) / abs(kkt)) < opts.tol_rel
        % Relative stopping criteria
        info.flag = 3;
        converged = true;
    elseif kkt < opts.tol_abs
        % Absolute stopping criteria
        info.flag = 4;
        converged = true;
    end
    old_kkt = kkt;
end

% Keep track of error history
if ~isempty(opts.errFcn)
    if isa(opts.errFcn,'function_handle')
        info.errHist(k+1) = opts.errFcn(x_k);
    elseif iscell(opts.errFcn)
        for i = 1:numel(opts.errFcn)
            fcn = opts.errFcn{i};
            info.errHist(k+1,i) = fcn(x_k);
        end
    end
end

if opts.storeIts 
    info.iterHist(k+1,:) = x_k;
end
% Fill up info if converged
if converged 
    if k == 0 
        info.flag = 9;
    end
    info.iter = k;
    info.f = f_k;
    %info.kkt = kkt;
    %info.msg = flag2msg(info.flag);
    if ~isempty(opts.errFcn)
        info.errHist = info.errHist(1:k+1,:);  
    end
    if opts.storeIts
        info.iterHist = info.iterHist(1:k+1,:);
    end
end
end % checkConvergence

function msg = flag2msg(flag)
assert(1 <= flag && flag <= 9)
% Just a list of human readable text strings to convert the flag return
% code into something readable by writing msg(flag) onto the screen.
msgs = {...
    'preprocessing';  % info.flag =  1
    'iterating';      % info.flag =  2
    'relative';       % info.flag =  3
    'absolute';       % info.flag =  4
    'stagnation';     % info.flag =  5
    'local minima';   % info.flag =  6
    'nondescent';     % info.flag =  7
    'maxlimit';       % info.flag =  8
    'x0 is sufficient'% info.flag =  9
};
msg = msgs{flag};
end % flag2msg