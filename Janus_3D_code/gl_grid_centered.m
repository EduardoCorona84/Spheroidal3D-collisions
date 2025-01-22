function [ u,v,rho ] = gl_grid_centered( p,r,xo,yo,zo)
%shifts GL grid
[u,v]=gl_grid(p);
rho=r*ones(size(u));
display(rho);

[x y z]=sph2cart(u,v,rho);


x=x+xo*ones(size(x));
y=y+yo*ones(size(y));
z=z+zo*ones(size(z));

[u,v,rho]=cart2sph(x,y,z);

end

