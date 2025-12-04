function [x, info] = callCVX(x0, A,b) %#ok<STOUT>
n = length(x0); %#ok<NASGU>
f = @(x) 1/2*dot(x,A*x) + dot(b,x); %#ok<NASGU>
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