function [ x, err,  iter, flag, convergence, msg] = ZERO_SR1( A, b, x0, max_iter, tol_rel, tol_abs, profile )
% Nic Rummel April 2025
if nargin<3
    x0 = zeros(N,1);
end
if nargin<4
    max_iter = floor(N/2);
end
if nargin<5
    tol_rel = 0.0001;
end
if nargin<6
    tol_abs = 10*eps; % Order of 10th of numerical precision seems okay
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
    'verbose',25, ...
    'nmax',max_iter, ...
    'tol',tol_rel, ...
    'x0', x0 ...
);

[x, iter, errStruct, ~] = zeroSR1(fcnGrad,[],h,prox,opts);
fK = errStruct(end,1);
fKm1 = errStruct(end-1,1);
gNorm =  errStruct(end,2);
err = min(abs(fK  - fKm1) / max([fK, fKm1, 1]), gNorm);

if iter == max_iter
    flag = 8;
    msg = 'maxlimit';
    return 
end

flag = 6;
msg = 'local minima';

convergence = [];
if profile
    % TODO, technically, we would want this match the above err...
    convergence = errStruct(:,2);
end

end

function hx = h_positiveOrthant(x)
if any(x<0)
    hx = Inf;
    return 
end
hx = 0;
end