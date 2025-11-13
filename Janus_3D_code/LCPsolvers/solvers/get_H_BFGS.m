function [H, opts] = get_H_BFGS(s, y, opts, bMask, debug)
n = numel(opts.b);
if ~exist('bMask', 'var') || isempty(bMask)
    bMask = false(n,1);
end
if ~exist('debug', 'var') || isempty(debug)
    debug = false;
end

[h0, rho, S, Y, opts] = updateQNMemory(s, y, opts);
B0 = @(x) x ./ h0;
H0 = @(x) x .* h0;
opts.prox.H0 = H0;
opts.prox.B0 = B0;

if isempty(S)
    H = H0; 
    opts.prox.U = [];
    opts.prox.V = [];
    opts.prox.B = B0; 
    return 
end


% Get a matrix free implementation of the inverse hessian approximation
H = @(g) apply_H(g, bMask, rho, S, Y, H0);
if nargout == 1
    return 
end
% If requested provide U and V such that B = B0 + U*U' + V*V'
% see N&W pg 184 for the unrolled formulas
% we don't use the compact representation because a sqrt may not exist of
% the inner matrix
r = size(S,2);
U = zeros(n, r);
V = zeros(n, r);
for i = 1:r
    s = S(~bMask,i); 
    y = Y(~bMask,i);
    U(~bMask,i) = y / sqrt(dot(s,y));
    v = B0(s);
    if i > 1 
        ucoeff = arrayfun(@(j) dot(U(~bMask,j),s), 1:i-1)';
        vcoeff = arrayfun(@(j) dot(V(~bMask,j),s), 1:i-1)';
        v = v + U(~bMask,1:i-1)*ucoeff - V(~bMask,1:i-1)*vcoeff;
    end
    V(~bMask,i) = v / sqrt(dot(s,v));
end
B = @(g) apply_B(g, bMask, rho, B0, U, V);
opts.prox.U = U; 
opts.prox.V = V;
opts.prox.B = B;
% Check the implementation
if debug 
    [Hexp, Bexp, Uexp, Vexp] = form_H_explicit(S(~bMask,:),Y(~bMask,:), ...
        H0, B0);
    Htest = eye(n);
    Btest = eye(n);
    for i = 1:n 
        Htest(:, i) = H(Htest(:,i));
        Btest(:, i) = B(Btest(:,i));
    end
    assert(norm(Htest(~bMask, ~bMask) - Hexp) / norm(Hexp) < 1e-6);
    assert(norm(Btest(~bMask, ~bMask)  - Bexp) / norm(Bexp) < 1e-6);
    assert(norm(U(~bMask,:) - Uexp) / norm(Uexp) < 1e-6);
    assert(norm(V(~bMask,:) - Vexp) / norm(Vexp) < 1e-6);
end 
end % get_H_BFGS

function p = apply_H(g, bMask, rho, S, Y, H0)
p = g(~bMask);
r = size(S,2);
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

p = H0(p);

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

function p = apply_B(g, bMask, rho, B0, U, V) %#ok<INUSD>

p = g(~bMask);
% r = size(S,2);
% alp = zeros(r,1);
% rs = zeros(r,1);
% for i = r:-1:1 
%     s = S(~bMask, i);
%     y = Y(~bMask, i);
%     if any(bMask)
%         rs(i) = 1/ dot(s,y);
%     else 
%         rs(i) = rho(i);
%     end
%     alp(i) = rs(i) * dot(s, p);
%     p = p - alp(i) * y;
% end

p = B0(p) + U*(U'*p) - V*(V'*p);

g(bMask) = 0;
g(~bMask) = p;
p = g;

end

function [H, B, U, V] = form_H_explicit(S,Y, H0, B0)
[n, r] = size(S);
U = zeros(n, r);
V = zeros(n, r);
H = H0(eye(n));
B = B0(eye(n));
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