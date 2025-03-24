function [f,g] = quadprog(x,Q,c,errFcn)
% f = qpnn(x,Q,c,errFcn)
%   returns the quadradic program objective problem 
%   to f(x) = 1/2*dot(x,Q*X) - dot(c,x)
% [f,g] = ...
%   return the gradient 
%   g(x) = Q*x-c
%
% [fHist,errHist] = normSquaredFunction()
%       will return the function history
%       (and error history as well, if errFcn was provided)
%       and reset the history to zero.
%   "fHist" is a record of f + extraFcn
%   (this is intended to be used where extraFcn is the non-smooth term "h")
%
% This function is (almost*) mathematically (not computationally) equivalent
%   to quadraticFunction( x, Q, c ) where
%   Q = A'*A and c = A'*b.
%   (*almost equivalent since there is a constant value difference in 
%    the objective function; you can use "constant" to change this)
%
% The Lipschitz constant of the gradient is 
%   the squared spectral norm of A, i.e., norm(A)^2
%
%
% March 4 2014, Stephen Becker, stephen.beckr@gmail.com
%
% See also quadraticFunction.m

persistent errHist fcnHist nCalls
if nargin == 0
   f = fcnHist(1:nCalls);
   g = errHist(1:nCalls);
   fcnHist = [];
   errHist = [];
   nCalls  = 0;
   return;
end

if nargin < 3 
    error('Must provide: x, Q,c')
end

if ~exist('errFcn','var') || isempty(errFcn)
    record = false;
elseif isa(errFcn, 'function_handle')
    record = true;
else 
    assert(isa(errFcn, 'bool'))
    record = errFcn;
end

f = 0.5*dot(x,Q*x) - dot(x,c);

if record
    nCalls = nCalls + 1;
    if length( errHist ) < nCalls
        % allocate more memory
        errHist(end:2*end) = 0;
        fcnHist(end:2*end) = 0;
    end
    fcnHist(nCalls) = f;
    if isa(errFcn, 'function_handle')
         errHist(nCalls) = errFcn(x);
    end
end

% ret grad
if nargout == 2
    g = Q*x - c;
end
