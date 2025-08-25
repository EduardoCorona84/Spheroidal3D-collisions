function [f, g, Ax] = quadraticLoss(x, A, b, Ax, Aq, eta)

if isempty(Ax) || isempty(Aq)
    Ax = A(x);
    f = 1/2*dot(x,Ax) + dot(b,x);
    if nargout == 1
        return 
    end
    g = Ax + b;
    return 
end
assert(exist('Aq','var') && exist('eta','var'))
assert(~isempty(Aq) && ~isempty(eta))
Ax = Ax + eta*Aq;
f = 1/2 *dot(x, Ax) + dot(x, b);
if nargout == 1
    return 
end
g = Ax + b; 

end