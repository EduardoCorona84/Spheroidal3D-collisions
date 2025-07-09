function [ x err  iter flag convergence msg] = APGD( A, b, x0, max_iter, tol_rel, tol_abs, profile )
% Copyright 2011, Kenny Erleben, DIKU

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

xk=x;
xkhat = ones(N,1);
yk=xk;
thetak=1;
Lk=norm(A*(xk-xkhat))/norm(xk-xkhat);
tk=1/Lk;
rmin=inf;

while (iter <= max_iter )
  Ayk=A*yk;
  g = Ayk + b;
  xkp1 = max(yk-tk*g,0);
  % Adjust the Lipshitz parameter
  tempvalue = 0.5*yk'*(Ayk)+yk'*b;
  Axkp1=A*xkp1;
  while (0.5*xkp1'*(Axkp1) + xkp1'*b >= tempvalue + g'*(xkp1-yk) + 0.5*Lk*norm(xkp1-yk))
      Lk=2*Lk;
      tk=1/Lk;
      xkp1 = max(yk-tk*g,0);
      Axkp1=A*xkp1;
  end
  
  % Nesterov 
  thetakp1 = 0.5*(-thetak*thetak + thetak*sqrt(thetak*thetak+4));
  betakp1=thetak*(1-thetak)/(thetak*thetak+thetakp1);
  ykp1=xkp1+betakp1*(xkp1-xk);
  
  r=rfallback(xkp1,Axkp1,b);
  if r<rmin
      rmin=r;
      xhat=xkp1;
  end
  
  % The APGD paper uses this as a convergence criteria, 
  % This is different from the merit function 0.5*(phi'*phi) used in BBPGD and mmNewton
  % To keep consistency, this criteria is not used
%    if r<tol_abs
%        flag = 4;
%        break;
%    end
  
  %--- Test all stopping criteria used ------------------------------------
  phi     = minmap(Axkp1 + b, xkp1);         % Calculate min map
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
  
  if g'*(xkp1-xk)>0
      ykp1 = xkp1;
      thetakp1=1;
  end
  Lk=0.9*Lk;
  tk=1/Lk;

  % Update iterate
  xk = xkp1;
  yk = ykp1;
  thetak=thetakp1;
  
  % Increment the number of iterations
  iter=iter+1;
end

if iter>=max_iter
  flag = 8;
  iter = iter - 1;
end

msg = msg{flag};
%{
disp('APGD result');
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

function [ r ] = rfallback(gamma,Ngamma,r)
gd=1e-6;
varphi = (gamma - max(gamma-gd*(Ngamma+r),0)) / (3*length(r)*gd) ;
r=norm(varphi);
end