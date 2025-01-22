function [ x err  iter flag convergence msg] = BBPGD( A, b, x0, max_iter, tol_rel, tol_abs, profile )
% May 2018 Wen Yan

% Just a list of human readable text strings to convert the flag return
% code into something readable by writing msg(flag) onto the screen.
msg = {'preprocessing';  % flag = 1
    'iterating';      % flag = 2
    'relative';       % flag = 3
    'absolute';       % flag = 4
    'stagnation';     % flag = 5
    'local minima';   % flag = 6
    'nondescent';     % flag = 7
    'maxlimit'        % flag = 8
    };

if nargin < 2
    error('Too few arguments');
end

N    = length(b); % Number of variables
flag = 1;

%--- Make sure we got good working default values -------------------------
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

%--- Make sure all values are valid ---------------------------------------
max_iter = max(max_iter,1);
tol_rel  = max(tol_rel,0);
tol_abs  = max(tol_abs,0);
x0       = max(0,x0);

%--- Here comes a bunch of magic constants --------------------------------

h       = 1e-7;    % Fixed constant used to evaluate the directional detivative
alpha   = 0.5;     % Step reduction parameter for projected Armijo backtracking line search
beta    = 0.001;   % Sufficent decrease parameter for projected Armijo backtracking line search
gamma   = 1e-28;   % Perturbation values used to fix near singular points in derivative
rho     = eps;     % Descent direction test parameter used to test if the Newton direction does a good enough job.

%--- Setup values need while iterating ------------------------------------

convergence = zeros(max_iter,1); % Used when profiling to measure the convergence rate

err     = Inf;         % Current error measure
x       = x0;          % Current iterate
iter    = 1;           % Iteration counter

flag    = 2;

% first step, plain GD
y = A*x+b;
stepsize = (y'*y) / (y'*(A*y+b));

while (iter <= max_iter )
    xold = x;
    yold = y;
    x = xold - stepsize*yold; % new x
    x(x<0)=0; % projection
    
    %--- Test all stopping criteria used ------------------------------------
    phi     = minmap(y,x);         % Calculate min map
    old_err = err;
    err     = 0.5*(phi'*phi);      % Natural merit function
    
    
    if profile
        convergence(iter) = err;
    end
    if (abs(err - old_err) / abs(old_err)) < tol_rel  % Relative stopping criteria
        flag = 3;
        break;
    end
    if err < tol_abs   % Absolute stopping criteria
        flag = 4;
        break;
    end
    
    % Test if the search direction is smaller than numerical precision. That is if it is too close to zero.
    if stepsize < 10*eps
        flag = 5;
        % Rather than just giving up we may just use the gradient direction
        % instead. However, I am lazy here!
        break;
    end
    
    y=A*x+b;
    
    % calculate B-B step size
    dx=x-xold;
    dy=y-yold;
    % BB1
    stepsize=(dx'*dx)/(dx'*dy);
    % BB2
    % stepsize=(dx'*dy)/(dy'*dy);
    % Increment the number of iterations
    iter=iter+1;
end

if iter>=max_iter
    flag = 8;
    iter = iter - 1;
end

msg = msg{flag};

%{
disp('BBPGD result');
disp(msg);
disp(iter);
disp(err);
%}
end

function [ phi ] = minmap(y,x)
% Auxiliary function used by the Minimum map Newton method
% Copyright 2011, Kenny Erleben, DIKU
phi  = min(y,x);
end
