% methods: 
%  - get_shc: pick easy spheroidal harmonic sigma, check that shc is the
%    same as shAna(sigma)
%  - get_sigma: pick identity vector for coefficients, confirm sigma
%    matches with shSyn(id)
%  - get_X: 
%      - Test translation: just shift spheroid
%      - Test rotation theta: just rotate by theta
%      - Test rotation phi: just rotate by phi
%  - get_X_targets
%  - separate_targets
%    plot

params=SpheroidalParameters;
params.a = 1;
params.u0 = [1.1 1.2 1.3];
params.centers= [0 0 0; 2 2 2; -1 -5 -1];
params.thetas = [0 pi/3 pi/8];
params.phis = [0 pi/4 6*pi/7];
params.p = 8;


% params.plot

[e,o]=test_X();
[e,o]=test_shc();
e
o
[e,o]=test_sigma();
e 
o

function [out, err] = test_shc()
    params = SpheroidalParameters;
    [th,ph] = gl_grid(8);
    f = Ynm(1,0,th,ph);
    params.sigma = f;
    fsh = shAna(f);
    params.get_shc;
    fsh_params = params.sigma_coefficients;
    err = log10(norm(fsh-fsh_params));
    out = err < -15;
end

function [out, err] = test_sigma()
    params = SpheroidalParameters;
    [th,ph] = gl_grid(8);
    f = Ynm(1,0,th,ph);
    
    fsh = zeros(9^2,1); fsh(3)=1; 

    params.sigma_coefficients = fsh;
    params.get_sigma;
    f_params = params.sigma;
    err = log10(norm(f-f_params));
    out = err < -15;
end

function [out, err] = test_X()
    % shift
    params = SpheroidalParameters;
    params.a = 1/2;
    params.u0 = [1.1 1.2];
    params.centers = [0 0 0; 1 1 1];
    params.phis=[0 0];
    params.thetas=[0 0];

    X1= prolate_spheroid_shape(8,params.u0(1),params.a,'cart');
    X2= prolate_spheroid_shape(8,params.u0(2),params.a,'cart');
    X2 = X2 + ones(size(X2));
    X = [X1; X2];

    getX=params.get_X;

    err = log10(norm(X-getX,2));
    out = err < -15;

    %rotate
    
end

