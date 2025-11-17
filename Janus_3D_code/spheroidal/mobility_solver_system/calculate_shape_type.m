function shape_type = calculate_shape_type(equ_radii, polar_radii)
    EQUALITY_TOL = 1e-14;
    ns = size(equ_radii, 1);
    shape_type = strings(ns, 1);
    for j=1:ns
        if abs(equ_radii - polar_radii) < EQUALITY_TOL
            shape_type(j) = 'sphere';
        elseif equ_radii > polar_radii
            shape_type(j) = 'oblate';
        else
            shape_type(j) = 'prolate';
        end
    end
end