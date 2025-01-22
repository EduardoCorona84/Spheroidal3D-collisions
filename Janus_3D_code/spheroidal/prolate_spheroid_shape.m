function S = prolate_spheroid_shape(p,u0,a,coords)
% Generate cartesian coordinates of order-p Gauss Legendre grid on a
% prolate spheroid surface fixed at u=u0, and scale 'a'.
if nargin ==3
    coords='cart';
end

[theta,phi]=gl_grid(p);

if strcmp(coords,'cart')
    x=a.*sqrt(u0^2-1).*sin(theta).*cos(phi);
    y=a.*sqrt(u0^2-1).*sin(theta).*sin(phi);
    z=a.*u0.*cos(theta);
    S=[x y z];
end

if strcmp(coords,'spheroidal')
S=[u0*ones(size(theta)) cos(theta) phi];

end

end