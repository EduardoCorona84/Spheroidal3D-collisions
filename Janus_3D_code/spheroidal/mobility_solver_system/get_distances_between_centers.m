function den = get_distances_between_centers(C)
    %{
    Calculates pairwise Euclidean distance between points in R^3.

    Inputs
    C - (double) n_b x 1 centers of bodies

    Outputs
    den - distances between bodies
    %}
    [Y_g1,  X_g1] = meshgrid(C(:,1), C(:,1));
    [Y_g2,  X_g2] = meshgrid(C(:,2), C(:,2));
    [Y_g3,  X_g3] = meshgrid(C(:,3), C(:,3));
    d1 = (X_g1 - Y_g1); d2 = (X_g2 - Y_g2); d3 = (X_g3 - Y_g3); 
    den = sqrt(d1.^2 + d2.^2 + d3.^2);

    % Diagonal entries are self-distances; only warn on zero off-diagonal distances.
    n = size(C, 1);
    if n > 1 && any(den(~eye(n)) == 0)
        warning('Zero distance between distinct centers in LOCAL_CenterDistance; handling of this is not implemented.');
    end
end
