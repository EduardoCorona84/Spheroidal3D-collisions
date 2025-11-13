function x = solveGeneralProx(y, B0) %#ok<STOUT>
n = length(y); %#ok<NASGU>
f = @(x) 1/2*dot(y-x, B0*(y-x)); %#ok<NASGU>
cvx_begin quiet
variable x(n)
minimize f(x)
subject to
0 <= x %#ok<NODEF,NOPRT>
cvx_end
end