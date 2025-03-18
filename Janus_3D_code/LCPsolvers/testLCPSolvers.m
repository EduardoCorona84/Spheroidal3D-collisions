%%
mfilePath = mfilename('fullpath');
if contains(mfilePath,'LiveEditorEvaluationHelper')
    mfilePath = matlab.desktop.editor.getActiveFilename;
end
[dirname, ~,~] = fileparts(mfilePath);
addpath(genpath(dirname))
%%
% build test problem so we can isolate just the solvers

% this problem is too hard right now
% N   = 12;
% A   = hilb(N);
% b   = ones(N,1);

% this problem is from the non-neg lsq on https://www.mathworks.com/help/matlab/ref/lsqnonneg.html
A = [0.0372    0.2869
     0.6861    0.7071
     0.6233    0.6245
     0.6344    0.6170];
 
b = [0.8587
     0.1781
     0.0747
     0.8405];
N   = size(A,2);

Q   = A'*A;
c   = A'*b;
normQ   = norm(Q);
%% use cvx to check that we have the correct answer
cvx_begin
        variable xRef(N)
        minimize sum_square(A*xRef-b)/2 
        subject to 
        0 <= xRef
cvx_end 
nrmXref = norm(xRef);
errFcn  = @(x) norm( x - xRef )/nrmXref;

% TODO: Note I think we can maybe make this slightly more efficient by not
% providing an upper bound...
prox          = @(x0,d,u,varargin) proj_rank1_box(0,Inf,x0,d,u);
h           = @(x) any(x < 0)*Inf; 
fcnGrad     = @(x) normSquaredFunction(x,A,[],b);

%% Solve with zeroSR1
% code from git@github.com:stephenbeckr/zeroSR1.git
opts = struct('N',N,'verbose',25,'nmax',4000,'tol',1e-13);
opts.L      = normQ; % optional
opts.errFcn = errFcn;

%  -- Default values usually fine --
% opts.BB     = true;
% opts.SR1_diagWeight=0.8;

tic
[xk,nit, errStruct,optsOut] = zeroSR1(fcnGrad,[],h,prox,opts);
% -- You can also call it this way, but can be slower --
% [xk,nit, errStruct,optsOut] = zeroSR1(fcnSimple,gradSimple,h,prox,opts);
tm = toc;
solverStr = 'zeroSR1';
fprintf('Final error for %15s is %.2e, took %.2f seconds\n', solverStr, errFcn(xk), tm );
figure(1); clf;
semilogy(errStruct(:,4) );
hold all
emphasizeRecent

%% Pure BB, no 0SR1
% code is also from the 0sr1 repo
opts.SR1 = false;
opts.BB_type    = 2;
opts.BB         = true;
[xk,nit, errStruct,optsOut] = zeroSR1(fcnGrad,[],h,prox,opts);
tm2 = toc;
solverStr = 'BB, no linesearch, i.e., basically SPG/SpaRSA';
fprintf('Final error for %15s is %.2e, took %.2f seconds\n', solverStr, errFcn(xk), tm2 );
semilogy(errStruct(:,4), '--' );

%% Compare to l-bfgs-b
% code from git@github.com:stephenbeckr/L-BFGS-B-C.git
opts    = struct( 'x0', zeros(N,1) );
opts.printEvery     = 400;
opts.m  = 5;
opts.errFcn     = @(x) errFcn(x);
% "outputFcn" will save values in the "info" output
opts.outputFcn  = opts.errFcn;
% Ask for very high accuracy
opts.pgtol      = 1e-14;
opts.factr      = 1e1;
opts.maxIts     = 5e4;
opts.maxTotalIts = 5e4;

f = @(x) 1/2 * norm(A*x - b)^2;
g = @(x) A'*(A*x-b); 
%driver_LeastSquares(x,A,b,offset);
tic
[x,fVal,info] = lbfgsb( {f,g} , zeros(N,1), Inf(N,1), opts );
toc
semilogy(info.err(:,3), '-.' );
%% PQN-L-BFGS
% code from https://www.cs.ubc.ca/~schmidtm/Software/PQN.html
options.order = -1;
options.optTol = 1e-10;
options.maxIter = 500;
options.corrections = 5;
options.bbInit = 0;
options.SPGoptTol = options.optTol;
options.SPGiters = options.maxIter;
options.maxProject = inf;
options.SPGtestOpt = 1;
options.errFcn = errFcn;
x0 = zeros(N,1);
function [ff,gg] = fg_tmp(x,A,b)
    ff = 1/2 * norm(A*x - b)^2;
    gg = A'*(A*x-b);
end
fg = @(x) fg_tmp(x,A,b);
funProj = @(x) max(0, x);
[xstar,fVal,funEvals, err] = minConF_PQN( fg ,x0,funProj,options);
semilogy(err(:,end), ':' );

% finish plotting
ylim([1e-12, 1e2])
ylabel('Error')
xlabel('Iteration')
title(['Nonnegative Least Squares Problem'])
hold all
emphasizeRecent
legend('zeroSR1','standard proximal gradient', 'l-bfgs-b', "pqn-l-bfgs");