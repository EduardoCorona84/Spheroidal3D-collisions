function diff = rel_kkt(x, A, b)
persistent olde
if ischar(x) && strcmpi(x,'reset')
    olde = [];
    diff = NaN;
    return
end
phi = min(x,A*x + b);
e = dot(phi, phi);
if isempty(olde)
    diff = NaN;
else
    diff = abs(e - olde) / abs(e);
end
olde = e;
end % rel_kkt