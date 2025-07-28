function [H, h0, U, V] = get_H_BFGS(k, s, y, opts, bMask, debug)

n = numel(s);
if ~exist('bMask', 'var') || isempty(bMask)
    bMask = false(n,1);
end

if ~exist('debug', 'var') || isempty(debug)
    debug = false;
end
n = numel(s);
[r, h0, rho, S, Y] = updateQNMemory(k, s, y, opts);
if r == 0
    H = @(g) g; 
    h0 = 1;
    U = [];
    V = [];
    return 
end

% Get a matrix free implementation of the inverse hessian approximation
H = @(g) apply_H(g, bMask, r, rho, S, Y, h0);

if nargout <= 2
    return 
end
% If requested provide U and V such that B = 1/h0 I + U*U' + V*V'
% see N&W pg 184 for the unrolled formulas
% we don't use the compact representation because a sqrt may not exist of
% the inner matrix
U = zeros(n, r);
V = zeros(n, r);
for i = 1:r
    s = S(~bMask,i); 
    y = Y(~bMask,i);
    U(~bMask,i) = y / sqrt(dot(s,y));
    v = s/h0 ;
    if i > 1 
        ucoeff = arrayfun(@(j) dot(U(~bMask,j),s), 1:i-1)';
        vcoeff = arrayfun(@(j) dot(V(~bMask,j),s), 1:i-1)';
        v = v + U(~bMask,1:i-1)*ucoeff - V(~bMask,1:i-1)*vcoeff;
    end
    V(~bMask,i) = v / sqrt(dot(s,v));
end

% Check the implmentation
if debug 
    [Hexp, Bexp, Uexp, Vexp] = form_H_explicit(S(~bMask,:),Y(~bMask,:), r, h0);
    Htest = eye(n);
    for i = 1:n 
        Htest(:, i) = H(Htest(:,i));
    end
    Btest = 1/h0*eye(n) + U*U' - V*V';
    assert(norm(Htest(~bMask, ~bMask) - Hexp) / norm(Hexp) < 1e-6);
    assert(norm(Btest(~bMask, ~bMask)  - Bexp) / norm(Bexp) < 1e-6);
    assert(norm(U(~bMask,:) - Uexp) / norm(Uexp) < 1e-6);
    assert(norm(V(~bMask,:) - Vexp) / norm(Vexp) < 1e-6);
end 
end % get_H_BFGS

function p = apply_H(g, bMask, r, rho, S, Y, h0)

p = g(~bMask);
alp = zeros(r,1);
rs = zeros(r,1);
for i = r:-1:1 
    s = S(~bMask, i);
    y = Y(~bMask, i);
    if any(bMask)
        rs(i) = 1/ dot(s,y);
    else 
        rs(i) = rho(i);
    end
    alp(i) = rs(i) * dot(s, p);
    p = p - alp(i) * y;
end

p = h0 * p;

for i = 1:r
    s = S(~bMask, i);
    y = Y(~bMask, i);
    b = rs(i) * dot(y,p);
    p = p + (alp(i) - b) * s;
end

g(bMask) = 0;
g(~bMask) = p;
p = g;

end

function [H, B, U, V] = form_H_explicit(S,Y, r, h0)
if ~exist('h0','var')
    h0 = 1;
end
[n, ~] = size(S);
U = zeros(n, r);
V = zeros(n, r);
H = eye(n) * h0;
B = eye(n) / h0;
for i = 1:r
    s = S(:,i); 
    y = Y(:,i); 
    rho = 1 / dot(y,s);
    if 1/rho < 1e-8
        % TODO something better than skip if the curvature condition is bad
        continue
    end
    UU = (eye(n) - rho*y*s');
    H = UU'*H*UU + rho*(s*s');

    v = B*s / sqrt(dot(s, B*s));
    u = y / sqrt(dot(y,s));

    B = B + u*u' - v*v' ;
    U(:, i) = u;
    V(:, i) = v;
end
% IDIOT CHECK 
assert(norm(B*H - eye(n)) < 1e-8)

end