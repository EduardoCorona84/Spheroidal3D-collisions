function [x, info] = callCVX(~, x0, opts) %#ok<STOUT>
n = length(x0); %#ok<NASGU>
% f = @(x) 1/2*sum_square(opts.A(x)+opts.b);
f = @(x) 1/2*dot(x,opts.AA*x) + dot(opts.b,x); %#ok<NASGU>
cvx_begin quiet
variable x(n)
minimize f(x)
subject to
0 <= x %#ok<NODEF,NOPRT>
cvx_end

info.iter = NaN;
info.kkt = NaN;
info.errHist = NaN;
info.iterHist = NaN;
end