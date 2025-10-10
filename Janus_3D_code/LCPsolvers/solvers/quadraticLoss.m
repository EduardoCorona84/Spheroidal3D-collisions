function [f, g] = quadraticLoss(x, A, b, Ax)
if ~exist('Ax', 'var') || isempty(Ax)
    Ax = A(x);
end

f = 1/2 *dot(x, Ax) + dot(x, b);

if nargout == 1
    return 
end
g = Ax + b; 

end