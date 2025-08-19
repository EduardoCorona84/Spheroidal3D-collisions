function surf_div = calculate_surf_div(params, S_self, sigma_x, sigma_y, sigma_z)
    %{
        Naively, one could expand out scalar densities in terms of
        spheroidal harmonics, and then see the action of the surface
        divergence operator on the expansion.

        Note that S_self is the vector of points for the associated
        spheroid in spheroidal coordinates.
    %}
    %%% Obtain SHC coefficients
    params.sigma = sigma_x; params.get_shc; shc_x = params.sigma_coefficients;
    params.sigma = sigma_y; params.get_shc; shc_y = params.sigma_coefficients;
    params.sigma = sigma_z; params.get_shc; shc_z = params.sigma_coefficients;

    p = params.p;
    sp = (p+1)^2;
    ii = (1:sp)'; nn=floor(sqrt(ii-1)); mm=ii-nn.^2-nn-1;

    %%% Set up coefficients for the surface divergence
    m = mm'; n = nn';
    u = S_self(:,1); v = S_self(:,2); phi = S_self(:,3);
    surf_div_coeffs = struct();
    if params.oblate
        surf_div_coeffs.shcX.Yn = (-1).*((-1).*(1+u.^2).*((-1)+v.^2)).^(-1/2).*(u.^2+v.^2).^(-2).*(((1+n).*v.^4+u.^4.*(2+((-1)+n).*v.^2)+u.^2.*(1+v.^2).*(1+n.*v.^2)).*cos(phi)+sqrt(-1).*m.*(u.^2+v.^2).^2.*sin(phi));
        surf_div_coeffs.shcX.Yn1 = ((-1)+m+(-1).*n).*v.*((-1)+v.^2).^(-1).*(((-1)+m+(-1).*n).^(-1).*(1+m+n).*(1+2.*n).*(3+2.*n).^(-1).*(1+u.^2).*((-1)+v.^2).*(u.^2+v.^2).^(-2)).^(1/2).*cos(phi);
        surf_div_coeffs.shcY.Yn = ((-1).*(1+u.^2).*((-1)+v.^2)).^(-1/2).*(u.^2+v.^2).^(-2).*(sqrt(-1).*m.*(u.^2+v.^2).^2.*cos(phi)+(-1).*((1+n).*v.^4+u.^4.*(2+((-1)+n).*v.^2)+u.^2.*(1+v.^2).*(1+n.*v.^2)).*sin(phi));
        surf_div_coeffs.shcY.Yn1 = ((-1)+m+(-1).*n).*v.*((-1)+v.^2).^(-1).*(((-1)+m+(-1).*n).^(-1).*(1+m+n).*(1+2.*n).*(3+2.*n).^(-1).*(1+u.^2).*((-1)+v.^2).*(u.^2+v.^2).^(-2)).^(1/2).*sin(phi);
        surf_div_coeffs.shcZ.Yn = u.*v.*(u.^2+v.^2).^(-2).*((-1)+((-1)+n).*u.^2+n.*v.^2);
        surf_div_coeffs.shcZ.Yn1 = ((-1)+m+(-1).*n).*((-1).*((-1)+m+(-1).*n).^(-1).*(1+m+n).*(1+2.*n).*(3+2.*n).^(-1)).^(1/2).*u.*(u.^2+v.^2).^(-1);
    else
        surf_div_coeffs.shcX.Yn = (u.^2+(-1).*v.^2).^(-2).*((-1).*((-1)+u.^2).*((-1)+v.^2)).^(-1/2).*(((-1).*(1+n).*v.^4+u.^4.*((-2)+v.^2+(-1).*n.*v.^2)+u.^2.*(1+v.^2).*(1+n.*v.^2)).*cos(phi)+(sqrt(-1)*(-1)).*m.*(u.^2+(-1).*v.^2).^2.*sin(phi));
        surf_div_coeffs.shcX.Yn1 = ((-1)+m+(-1).*n).*v.*(((-1)+m+(-1).*n).^(-1).*(1+m+n).*(1+2.*n).*(3+2.*n).^(-1).*((-1)+u.^2).*((-1)+v.^2).^(-1)).^(1/2).*((-1).*u.^2+v.^2).^(-1).*cos(phi);
        surf_div_coeffs.shcY.Yn = (-1).*(u.^2+(-1).*v.^2).^(-2).*((-1).*((-1)+u.^2).*((-1)+v.^2)).^(-1/2).*((sqrt(-1)*(-1)).*m.*(u.^2+(-1).*v.^2).^2.*cos(phi)+((1+n).*v.^4+u.^4.*(2+((-1)+n).*v.^2)+(-1).*u.^2.*(1+v.^2).*(1+n.*v.^2)).*sin(phi));
        surf_div_coeffs.shcY.Yn1 = ((-1)+m+(-1).*n).*v.*(((-1)+m+(-1).*n).^(-1).*(1+m+n).*(1+2.*n).*(3+2.*n).^(-1).*((-1)+u.^2).*((-1)+v.^2).^(-1)).^(1/2).*((-1).*u.^2+v.^2).^(-1).*sin(phi);
        surf_div_coeffs.shcZ.Yn = u.*(((-1)+v.^2).*((-1).*u.^2+v.^2).^(-1)).^(1/2).*((-1).*u.^2.*v.*(1+(-1).*v.^2).^(-1/2).*(u.^2+(-1).*v.^2).^(-3/2)+v.*(u.^2+(-1).*v.^2).^(-2).*(((-1)+v.^2).^(-1).*((-1).*u.^2+v.^2)).^(1/2)+n.*v.*(((-1)+v.^2).*((-1).*u.^2+v.^2)).^(-1/2));
        surf_div_coeffs.shcZ.Yn1 = ((-1)+m+(-1).*n).*u.*((-1).*((-1)+m+(-1).*n).^(-1).*(1+m+n).*(1+2.*n).*(3+2.*n).^(-1).*(u.^2+(-1).*v.^2).^(-2)).^(1/2);
    end

    %%% Get spheroidal harmonics
    nt = length(u);
    Yr=zeros(2*nt,sp);

    for n=0:p
        %%% Store Yn harmonics
        Yn=Ynm(n, [], real(acos(v)), phi);
        Yr(1:nt, n^2+1:(n+1)^2)=Yn;

        %%% Store Yn+1 harmonics
        Yn1=Ynm(n+1, -n:n, real(acos(v)), phi);
        Yr(nt+1:2*nt, n^2+1:(n+1)^2) = Yn1;
    end

    %%% Finally, compute the final surface divergence term.
    common_coeff = 1 / params.a;
    shc_X_coeffs = surf_div_coeffs.shcX.Yn.*Yr(1:nt,:) + surf_div_coeffs.shcX.Yn1.*Yr(nt+1:end,:);
    shc_Y_coeffs = surf_div_coeffs.shcY.Yn.*Yr(1:nt,:) + surf_div_coeffs.shcY.Yn1.*Yr(nt+1:end,:);
    shc_Z_coeffs = surf_div_coeffs.shcZ.Yn.*Yr(1:nt,:) + surf_div_coeffs.shcZ.Yn1.*Yr(nt+1:end,:);
    surf_div = (common_coeff.' .* shc_X_coeffs) * shc_x + ...
               (common_coeff.' .* shc_Y_coeffs) * shc_y + ...
               (common_coeff.' .* shc_Z_coeffs) * shc_z;
end