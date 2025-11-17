function [u0, a] = calculate_u0_and_a_from_radii(shape_type, equ_radii, polar_radii)
%{
Helper function to determine the parameters u0 and a that are associated with spheroids.
%}
    if isempty(shape_type)
        shape_type = calculate_shape_type(equ_radii, polar_radii)
    end

    if strcmp(shape_type, 'prolate')
        a = sqrt(polar_radii.^2 - equ_radii.^2);
        u0 = polar_radii ./ a;
    elseif strcmp(shape_type, 'oblate')
        a = sqrt(equ_radii.^2 - polar_radii.^2);
        u0 = equ_radii ./ a;
    else
        error('Invalid shape passed in: should be prolate or oblate.');
    end
end