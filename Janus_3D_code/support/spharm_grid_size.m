function [nu,nv, p]=spharm_grid_size(p,ntot)
% SPHARM_GRID_SIZE returns the lattitude and longitude grid size for
% spherical harmonic order p.

if(nargin<2), ntot=-1;end
if(isempty(p))
    p  = (round(sqrt(2*ntot+1))-1)/2;
end
nu = p + 1;
nv = 2*p;

if(nargin>1), assert(nv*nu==ntot);end