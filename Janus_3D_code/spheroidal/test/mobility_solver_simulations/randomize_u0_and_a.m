%{
Utility code to pick out u0's and a's for simulation code.

The reasoning for choosing a = 1./u0 is to force the discretization points
to be somewhat close to each other so we don't have to up the harmonic order
in order to get the same level of accuracy.
%}

function [u0, a, oblate, centers] = randomize_u0_and_a(ns, max_u0, min_u0, randomize_oblate_flag)
    %{
    Inputs
        ns : strictly positive integer
            number of bodies
        max_u0 : float
        min_u0 : float
        randomize_oblate_flag : boolean
            by default everything is prolate; set this to be true if one wants a mix of oblates/prolates

    Outputs
        u0 : 1 x ns numeric array
        a : 1 x ns numeric array
        oblate : 1 x ns logical array
            array to denote whether something is an oblate or not
    %}
    assert(ns > 0, 'Number of bodies must be strictly positive.');
    u0 = zeros(1, ns);
    a = zeros(1, ns);

    %% u0, a, oblate randomization
    if randomize_oblate_flag
        oblate = randi([0 1], [1 ns]);
    else
        oblate = false(1, ns);
    end

    for i=1:ns
        u0(i) = rand(1, ns)*max_u0 + min_u0;
        if oblate(i)
            a(i) = 1/u0(i);
        else
            a(i) = 1/sqrt(1 + u0(i)^2);
        end
    end

    %% Centers randomization
    % The idea is to ensure everything is a major-radius away
    % In this case, the major radii of all spheroid is 1.
end