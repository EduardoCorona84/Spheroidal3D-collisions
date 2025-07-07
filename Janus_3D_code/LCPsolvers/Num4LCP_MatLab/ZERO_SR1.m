function [ x, err,  iter, flag, convergence, msg] = ZERO_SR1( A, b, x0, max_iter, tol_rel, tol_abs, profile )
% Nic Rummel April 2025
if nargin<3
    x0 = zeros(N,1);
end
if nargin<4
    max_iter = floor(N/2);
end
if nargin<5
    tol_rel = 1e-6;
end
if nargin<6
    tol_abs = 1e-4; 
end
if nargin<7
    profile = true;
end
% set up this solvers options
prox    = @(x,d,u,varargin) proj_rank1_box(0, Inf, x,d,u); 
h       = @(x)h_positiveOrthant(x);
fcnGrad = @(x) quadprog(x,A,-b);
opts = struct( ...
    'N',size(A,2), ...
    'verbose',false, ...
    'nmax',max_iter, ...
    'tol', 0, ...
    'tol_rel',tol_rel, ...
    'tol_abs',tol_abs, ...
    'x0', x0, ...
    'errFcn', @(x) NaN, ...
    'Q', A ...
);

[x, iter, errStruct, ~] = zeroSR1(fcnGrad,[],h,prox,opts);
fK = errStruct.f(end);
fKm1 = errStruct.f(end-1);
gNorm =  errStruct.gnorm(end);
err = min(abs(fK  - fKm1) / max([fK, fKm1, 1]), gNorm);

convergence = [];
if profile
    % TODO, technically, we would want this match the above err...
    convergence = errStruct;
end

if iter == max_iter
    flag = 8;
    msg = 'maxlimit';
    return 
end

flag = 6;
msg = 'local minima';

end

function hx = h_positiveOrthant(x)
if any(x<0)
    hx = Inf;
    return 
end
hx = 0;
end