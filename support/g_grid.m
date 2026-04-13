function [qp,qw] = g_grid(n)
% G_GRID(n) - Returns the n-point Gauss-Legendre quadrature nodes and weights
%
% qp is the vector of quadrature points in [-1,1] interval and qw is the
% quadrature weight vector.
%
% This code is courtesy of L.N. Trefethen, Spectral method in Matlab, page 129.
%
% SEE ALSO: GL_GRID  
  
if(nargin==0), testThis();return;end

beta = .5./sqrt(1-(2*(1:n-1)).^(-2));
T = diag(beta,1) + diag(beta,-1);
[V,D] = eig(T);
qp = diag(D);
[qp,ind] = sort(qp); qp = qp';
qw = 2*V(1,ind).^2;

function testThis()
% TESTTHIS - Tester function for grule
%   
  m = 12;
  [qp,qw] = g_grid(m);
  fprintf(' With %d points,\n polynomial of order of...\n',m);
  fprintf(' -----------------------\n');
  for ii=[1:4:2*m-2 2*m-1:2*m+3]
    p = rand(1,ii+1);
    intN = sum(qw.*polyval(p,qp));
    
    p = [p./(ii+1:-1:1) 0];
    intA = polyval(p,1)-polyval(p,-1);
    fprintf(' %-2g has error: %-2.3e\n',ii,abs(1-intN/intA));
	if(ii==2*m-1), fprintf(' -----------------------\n');end
  end
  fprintf(' -----------------------\n');