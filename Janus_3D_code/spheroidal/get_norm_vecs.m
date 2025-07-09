%{
    Extracts the normal vector given p and u0 for a spheroid.

    Note that this is from SpheroidalParameters.m
%}


function nvecs = get_norm_vecs(p, u0, oblate_flag)
    np = 2*p*(p+1);
    [theta_x, phi_x] = gl_grid(p);
    v_x=cos(theta_x);
    if ~oblate_flag
        nvecs = 1./sqrt((u0.^2-1).*(u0.^2-v_x.^2)).*[u0.*sqrt(u0.^2-1).*sqrt(1-v_x.^2).*cos(phi_x),...
            u0.*sqrt(u0.^2-1).*sqrt(1-v_x.^2).*sin(phi_x),...
            (u0.^2-1).*v_x];
    else
        nvecs = 1./sqrt((u0.^2+1).*(u0.^2+v_x.^2)).*[u0.*sqrt(u0.^2+1).*sqrt(1-v_x.^2).*cos(phi_x),...
            u0.*sqrt(u0.^2+1).*sqrt(1-v_x.^2).*sin(phi_x),...
            (u0.^2+1).*v_x];
    end
end