%%
rng(1)
n = 10;
m = 2;
A = randn(n,n) + 10*eye(n);
A = A+A';
b = randn(n,1);
B = eye(m,n);%rand(m,n);
c = randn(m,1);
while all(B * (A \ b) > c) || all(B * (A \ b) <= c)
    c = randn(m,1);
end
%% Unconstrained Problem
disp('===========================')
disp('== Unconstrained Problem ==')
% Primal is just qp
cvx_begin quiet
variable x(n)
minimize 1/2*dot(x,A*x) + dot(b,x)
cvx_end
disp('Primal Rel Err compared to linSlv')
disp(norm(x + A\b) / norm(A\b))
% Dual is also qp but with inverse
cvx_begin quiet
variable u(n)
minimize 1/2*dot(-u-b,A\(-u-b))
subject to 
u == 0
cvx_end
xViaDual =  A \ (-u - b);
disp('Primal Rel Err Lagrange')
disp(norm(xViaDual - x) / norm(x))
%% Constrained Problem 
disp('===========================')
disp('=== Constrained Problem ===')
%% Primal 
cvx_begin quiet
variable x(n)
dual variable u
minimize 1/2*dot(x,A*x) + dot(b,x)
subject to 
    u : c <= B *x  %#ok<NOPTS>
cvx_end
%% Dual Lagrange
cvx_begin quiet
variable uL(m)
dual variable xF
minimize 1/2* dot(A \ B'*uL , B'*uL) -  dot(uL, B*(A\b) + c)
subject to 
    xF : 0 <= uL
cvx_end
xL =  A \ (B' * uL - b);
disp('Primal Rel Err Lagrange')
disp(norm(xL - x) / norm(x))
disp('Dual Rel Err Lagrange')
disp(norm(uL - u) / norm(u))
%% Dual Fenchel
cvx_begin quiet
variable uF(m)
dual variable xF
minimize 1/2* dot(A \ (B'*uF + b) , B'*uF + b)
subject to 
    xF : uF <= c
cvx_end
xF =   - A \ (B' * uF + b);
disp('Primal Rel Err 1')
disp(norm(xF - x) / norm(x))
disp('Dual Rel Err 1')
disp(norm(uF - u) / norm(u))