function [u,v] = gl_grid(p)
% GL_GRID - Returns the Gauss-Legendre--uniform grid on the unit sphere
% [0,pi]x[0,2pi).
%
% SEE ALSO: G_GRID
%
  
[nu,nv]  = spharm_grid_size(p);

lambda = (0:nv-1)'*2*pi/nv;
theta  = acos(g_grid(nu));
[v, u] = meshgrid(lambda,theta);
u = u(:); v = v(:);