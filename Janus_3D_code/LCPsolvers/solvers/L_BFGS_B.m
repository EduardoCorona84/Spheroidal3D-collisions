function [ x, info] = L_BFGS_B(fg, x0, opts )
% Nic Rummel April 2025
opts = defaultLCPOpts(opts, x0);


% The error bounds are slightly different for this code than that of the other
% solvers in particular
%  - 'fctr' parameter specifies to break when
%    |f_{k+1} - f_{k}| / max(f_{k+1}, f_{k}, 1) < factr*eps()
%    so we choose to make this as loose as possible via tol_rel & tol_abs
%  - 'pgtol' parameter specifies to break when
%    max{|proj g_i | i = 1, ..., n} <= pgtol
%    this is effectively a bound on the inf-norm of the gradient, so we use
%    the abs_tol
fctr = 1; %max(opts.tol_rel, opts.tol_abs) / eps();
lbfgsOpts = struct(...
    'x0', x0, ...
    'printEvery', Inf, ...
    'm', opts.qn.m, ...
    'pgtol', opts.kkt_abs, ...
    'factr', fctr, ...
    'maxIts', opts.max_iter, ...
    'maxTotalIts', opts.linesearch.budget*opts.max_iter);
lbfgsOpts.errFcn = opts.errFcn;
if ishandle(opts.errFcn)
        info.errHist = opts.errFcn(x0);
elseif iscell(opts.errFcn)
    info.errHist = zeros(1,numel(opts.errFcn)+1);
    for i = 1:numel(opts.errFcn)
        fcn = opts.errFcn{i};
        info.errHist(1,i) = fcn(x0);
    end
    info.errHist(1,end) = fg(x0, []);
end
n = numel(x0);
lb = zeros(n,1);
ub = Inf(n,1);
this_fg = @(x) fg(x, []);
[x, ~, lbfgsInfo] = lbfgsb(this_fg , lb, ub, lbfgsOpts );
info.iter = lbfgsInfo.iterations;
if ~isempty(opts.errFcn)
    info.errHist = cat(1, info.errHist, cat(2, lbfgsInfo.err(:,3:end), lbfgsInfo.err(:,1)));
    info.kkt = lbfgsInfo.err(end,3); % assumes that the first errFcn is kkt
else 
    [~,g] = this_fg(x);
    phi = min(x,g);
    info.kkt = 1/2 * dot(phi,phi);
end

if info.iter == opts.max_iter
    info.flag = 8;
    info.msg = 'maxlimit';
    return 
end

info.flag = 6;
info.msg = 'local minima';
info.iterHist = NaN;

end