function f=Yp(n,m,p)
%spheroidal coordinates of the grid points
[u,v]=gl_grid(p);
f=Ynm(n,m,u,v);
end
