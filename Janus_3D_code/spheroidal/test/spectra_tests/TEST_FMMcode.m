params_source = SpheroidalParameters;
params_source.isReal = true;
params_source.u0 = 2/sqrt(3); % AR = 2
params_source.a = 1/params_source.u0;
params_source.oblate = false;
params_source.matvec_eta = 1;

p = 16;
np = 2*p*(p+1);

[u, v] = gl_grid(p);
params_source.sigma = cos(u) .* cos(v);

params_target = SpheroidalParameters;
params_target.centers = [5 0 0];
params_target.isReal = true;
params_target.u0 = 3/sqrt(5); % AR = 1.5
params_target.a = 1/params_target.u0;
params_target.oblate = false;
params_target.sigma = params_source.sigma;

X_src = params_source.get_X();
X_trg = params_target.get_X();
nu = params_target.get_Norm();

matvec_res_with_FMM = spheroidalMatVec(params_source, 'SL', X_trg, nu, true);
matvec_res_with_smoothquad = spheroidalMatVec(params_source, 'SL', X_trg, nu, false);

abs_err = norm(matvec_res_with_smoothquad - matvec_res_with_FMM)

matvec_res_with_FMM = spheroidalMatVec(params_source, 'DL', X_trg, nu, true);
matvec_res_with_smoothquad = spheroidalMatVec(params_source, 'DL', X_trg, nu, false);

abs_err = norm(matvec_res_with_smoothquad - matvec_res_with_FMM)

matvec_res_with_FMM = spheroidalMatVec(params_source, 'SP', X_trg, nu, true);
matvec_res_with_smoothquad = spheroidalMatVec(params_source, 'SP', X_trg, nu, false);

abs_err = norm(matvec_res_with_smoothquad - matvec_res_with_FMM)

matvec_res_with_FMM = spheroidalMatVec(params_source, 'DP', X_trg, nu, true);
matvec_res_with_smoothquad = spheroidalMatVec(params_source, 'DP', X_trg, nu, false);

abs_err = norm(matvec_res_with_smoothquad - matvec_res_with_FMM)