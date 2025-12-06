function [converged, info] = checkConvergence(k, f_k, x_k, ...
    grad_k, eta, info, opts, x_km1_, grad_km1_)
persistent old_kkt x_km1;
if k == 0
    old_kkt = [];
    x_km1 = [];
else
    if exist('x_km1_','var') && ~isempty(x_km1_) 
        x_km1 = x_km1_;
    end
    if exist('grad_km1_','var') && ~isempty(grad_km1_)
        phi = min(grad_km1_,x_km1);
        old_kkt = dot(phi, phi);
    end
end
% Check for convergence
converged = false;
%% kkt conditions / LCP being satisified is equivalent to
% the grad_k(i) = 0 \perp x_k(i) = 0 (for symetric A)
phi = min(grad_k,x_k);
kkt = dot(phi, phi);
%% If the sequence is a cauchy sequence we are converging.
iterErr = Inf;
if k > 0
    iterErr = norm(x_k - x_km1);
end
%% 
if k >= opts.max_iter
    info.flag =  8;
    converged = true;
elseif ~isempty(opts.kkt_rel) && ~isempty(old_kkt) && (abs(kkt - old_kkt) / abs(kkt)) < opts.kkt_rel
    % Relative stopping criteria
    info.flag = 3;
    converged = true;
elseif ~isempty(opts.kkt_abs) && kkt < opts.kkt_abs
    % Absolute stopping criteria
    info.flag = 4;
    converged = true;
elseif ~isempty(opts.arg_rel) && iterErr / max(norm(x_k), norm(x_km1)) < opts.arg_rel
    % Relative stopping criteria on the iterates
    info.flag = 5;
    converged = true;
elseif ~isempty(opts.arg_abs) && iterErr < opts.arg_abs
    % Absolute stopping criteria on the iterates
    info.flag = 6;
    converged = true;
elseif ~isempty(opts.step_abs) && eta < 10*eps
    % if step direction gets too small give up
    info.flag =  7;
    converged = true;
end
old_kkt = kkt;
x_km1 = x_k;
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
    info.kkt = kkt;
    info.msg = flag2msg(info.flag);
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
    'relative kkt';   % info.flag =  3
    'absolute kkt';   % info.flag =  4
    'relative arg';   % info.flag =  5
    'absolute arg';   % info.flag =  6
    'stagnation';     % info.flag =  7
    'maxlimit';       % info.flag =  8
    'x0 is sufficient'% info.flag =  9
};
msg = msgs{flag};
end % flag2msg