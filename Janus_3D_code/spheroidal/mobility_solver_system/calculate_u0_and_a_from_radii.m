function [u0, a] = calculate_u0_and_a_from_radii(shape_type, equatorial_radius, polar_radius)
%{
Helper function to determine the parameters u0 and a that are associated with spheroids.
%}
    if strcmp(shape_type, 'prolate')
        a = sqrt(polar_radius.^2 - equatorial_radius.^2);
        u0 = polar_radius ./ a;
    elseif strcmp(shape_type, 'oblate')
        a = sqrt(equatorial_radius.^2 - polar_radius.^2);
        u0 = equatorial_radius ./ a;
    else
        error('Invalid shape passed in: should be prolate or oblate.');
    end
end