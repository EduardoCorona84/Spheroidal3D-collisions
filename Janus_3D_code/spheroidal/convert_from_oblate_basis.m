function [X, Y, Z] = convert_from_oblate_basis(U_comp, V_comp, PHI_comp, u, v, phi)
    %{
        Given a vector in the spheroidal basis; that is,
            v = A\hat{u} + B\hat{v} + C\hat{\phi},
        convert the vector to Cartesian coordinates so that
            v = X\hat{x} + Y\hat{y} + Z\hat{z},
        where we use the canonical Cartesian basis.
    %}
    u2 = u.^2;
    v2 = v.^2;

    sqrt_u2m1 = sqrt(u2 + 1);
    sqrt_1mv2 = sqrt(1 - v2);
    sqrt_u2mv2 = sqrt(u2 + v2);
    
    cos_phi = cos(phi);
    sin_phi = sin(phi);
    
    M11 = (u .* sqrt_1mv2 ./ sqrt_u2mv2) .* cos_phi;
    M21 = (u .* sqrt_1mv2 ./ sqrt_u2mv2) .* sin_phi;
    M31 = (v .* sqrt_u2m1 ./ sqrt_u2mv2);
    
    M12 = (-v .* sqrt_u2m1 ./ sqrt_u2mv2) .* cos_phi;
    M22 = (-v .* sqrt_u2m1 ./ sqrt_u2mv2) .* sin_phi;
    M32 = (u .* sqrt_1mv2 ./ sqrt_u2mv2);
    
    M13 = -sin_phi;
    M23 =  cos_phi;
    M33 = zeros(size(u));
    
    X = M11.*U_comp + M12.*V_comp + M13.*PHI_comp;
    Y = M21.*U_comp + M22.*V_comp + M23.*PHI_comp;
    Z = M31.*U_comp + M32.*V_comp + M33.*PHI_comp;
end