function shape_type = calculate_shape_type(equ_radii, polar_radii)
    EQUALITY_TOL = 1e-14;
    equ_radii = equ_radii(:);
    polar_radii = polar_radii(:);
    assert(numel(equ_radii) == numel(polar_radii), ...
        'equ_radii and polar_radii must have the same number of elements.');

    ns = numel(equ_radii);
    shape_type = strings(ns, 1);
    for j=1:ns
        if abs(equ_radii(j) - polar_radii(j)) < EQUALITY_TOL
            shape_type(j) = 'sphere';
        elseif equ_radii(j) > polar_radii(j)
            shape_type(j) = 'oblate';
        else
            shape_type(j) = 'prolate';
        end
    end
end
