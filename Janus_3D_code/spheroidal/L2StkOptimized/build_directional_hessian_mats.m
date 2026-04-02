function [Umat, Vmat, PHImat] = build_directional_hessian_mats(region, nu, a, oblate)
%{
Build the precombined U/V/PHI blocks for the directional Hessian.

Usage:
  [Umat, Vmat, PHImat] = build_directional_hessian_mats(region, nu_cart, a, oblate)

Inputs
  region : region struct from spheroidal_cartesian_target_precompute
  nu     : nt x 3 Cartesian normal components
  a      : focal length
  oblate : oblate/prolate geometry flag

Notes
  - The outputs are already precombined against cached data, so the
    apply stage is just
      HnU   = Umat   * Gshc
      HnV   = Vmat   * Gshc
      HnPHI = PHImat * Gshc
%}
nu_sph = cartNu2spheroidal(nu, [region.u, region.v, region.phi], a, oblate);

sp = size(region.Yr0, 2);
ii = (1:sp)';
nn = floor(sqrt(ii - 1));
mm = ii - nn.^2 - nn - 1;

nrow = (nn + 1).';
n2row = (nn + 2).';
mrow = mm.';
m2row = (mm.^2).';
ones_row = ones(1, sp);

alpha_row = ((mm - nn - 1) .* sqrt(((nn + mm + 1) .* (2 * nn + 1)) ./ ...
    ((nn - mm + 1) .* (2 * nn + 3)))).';
beta_row = ((mm - nn - 2) .* sqrt(((nn + mm + 2) .* (2 * nn + 3)) ./ ...
    ((nn - mm + 2) .* (2 * nn + 5)))).';

u = region.u(:);
v = region.v(:);
nu_u = nu_sph(:, 1);
nu_v = nu_sph(:, 2);
nu_phi = nu_sph(:, 3);

if oblate
    Delta = u.^2 + v.^2;
    root_u = sqrt(u.^2 + 1);
    denom_u = u.^2 + 1;
else
    Delta = u.^2 - v.^2;
    root_u = sqrt(u.^2 - 1);
    denom_u = u.^2 - 1;
end
root_v = sqrt(1 - v.^2);
root_Delta = sqrt(Delta);
denom_v = 1 - v.^2;

common_coeffs = region.common_coeffs.';
Fr = region.Fr;
Fp = region.Fp;
Fpp = region.Fpp;
Yr0 = region.Yr0;
Yr1 = region.Yr1;
Yr2 = region.Yr2;

if oblate
    %% Umat
    U_coef0_f = (1i .* nu_phi .* u ./ (denom_u .* root_Delta .* root_v)) * mrow ...
                + (nu_v .* u .* v .* root_u ./ (Delta.^2 .* root_v) - nu_u .* v.^2 ./ Delta.^2) * nrow;

    U_coef0_fp = (nu_phi ./ (root_Delta .* root_v)) * mrow ...
                + (1i .* nu_v .* v .* root_u .* root_v ./ Delta.^2 + 1i .* nu_u .* u .* denom_v ./ Delta.^2) * ones_row ...
                - (1i .* nu_v .* root_u .* v ./ (Delta .* root_v)) * nrow;

    U_coef0_fpp = (nu_u .* denom_u ./ Delta) * ones_row;

    U_coef1_f = (nu_v .* u .* root_u ./ (Delta.^2 .* root_v) - nu_u .* v ./ Delta.^2) * alpha_row;
    U_coef1_fp = (-1i .* nu_v .* root_u ./ (Delta .* root_v)) * alpha_row;

    Umat = -common_coeffs .* ( ...
        (U_coef0_f .* Fr + U_coef0_fp .* Fp + U_coef0_fpp .* Fpp) .* Yr0 + ...
        (U_coef1_f .* Fr + U_coef1_fp .* Fp) .* Yr1);

    %% Vmat
    V_coef0_f = (-1i .* nu_phi .* v ./ (root_u .* root_Delta .* denom_v)) * (mrow .* n2row) ...
                - (nu_v ./ (Delta.^2 .* denom_v)) .* (nrow .* (u.^2 .* (1 + nrow .* v.^2) + n2row .* v.^4)) ...
                + (nu_u .* root_u .* u .* v ./ (Delta.^2 .* root_v)) * nrow;

    V_coef0_fp = (-1i .* nu_v .* u .* denom_u ./ Delta.^2) * ones_row ...
                + (1i .* nu_u .* root_u .* v ./ (Delta.^2 .* root_v)) .* (denom_v - nrow .* Delta);

    V_coef1_f = (-1i .* nu_phi ./ (root_u .* root_Delta .* denom_v)) * (mrow .* alpha_row) ...
                - (nu_v .* v ./ (Delta.^2 .* denom_v)) .* alpha_row .* (-1 + 2 .* (nrow + 1) .* u.^2 + (2 .* nrow + 3) .* v.^2) ...
                + (nu_u .* root_u .* u ./ (Delta.^2 .* root_v)) * alpha_row;

    V_coef1_fp = (-1i .* nu_u .* root_u ./ (Delta .* root_v)) * alpha_row;
    V_coef2_f = -(nu_v ./ (Delta .* denom_v)) * (alpha_row .* beta_row);

    Vmat = -common_coeffs .* ( ...
        (V_coef0_f .* Fr + V_coef0_fp .* Fp) .* Yr0 + ...
        (V_coef1_f .* Fr + V_coef1_fp .* Fp) .* Yr1 + ...
        (V_coef2_f .* Fr) .* Yr2);

    %% PHImat
    PHI_coef0_fp = (-1i .* nu_phi .* u ./ Delta) * ones_row ...
                    + (nu_u ./ (root_v .* root_Delta)) * mrow;

    PHI_coef0_f = (nu_phi ./ (denom_u .* denom_v)) * m2row ...
                + (nu_phi .* v.^2 ./ (denom_v .* Delta)) * nrow ...
                - (1i .* nu_v .* v ./ (root_u .* denom_v .* root_Delta)) * (mrow .* n2row) ...
                + (1i .* nu_u .* u ./ (denom_u .* root_v .* root_Delta)) * mrow;

    PHI_coef1_f = (nu_phi .* v ./ (denom_v .* Delta)) * alpha_row ...
                - (1i .* nu_v ./ (root_u .* denom_v .* root_Delta)) * (mrow .* alpha_row);

    PHImat = -common_coeffs .* ( ...
        (PHI_coef0_f .* Fr + PHI_coef0_fp .* Fp) .* Yr0 + ...
        (PHI_coef1_f .* Fr) .* Yr1);
else
    %% Umat
    U_coef0_f = (-1i .* nu_phi .* u ./ (denom_u .* root_Delta .* root_v)) * mrow ...
                - (nu_v .* u .* v .* root_u ./ (Delta.^2 .* root_v) + nu_u .* v.^2 ./ Delta.^2) * nrow;

    U_coef0_fp = (1i .* nu_phi ./ (root_Delta .* root_v)) * mrow ...
                + (nu_v .* v .* root_u .* root_v ./ Delta.^2 + nu_u .* u .* denom_v ./ Delta.^2) * ones_row ...
                + (nu_v .* root_u .* v ./ (Delta .* root_v)) * nrow;

    U_coef0_fpp = (nu_u .* denom_u ./ Delta) * ones_row;

    U_coef1_f = -(nu_v .* u .* root_u ./ (Delta.^2 .* root_v) + nu_u .* v ./ Delta.^2) * alpha_row;
    U_coef1_fp = (nu_v .* root_u ./ (Delta .* root_v)) * alpha_row;

    Umat = common_coeffs .* ( ...
        (U_coef0_f .* Fr + U_coef0_fp .* Fp + U_coef0_fpp .* Fpp) .* Yr0 + ...
        (U_coef1_f .* Fr + U_coef1_fp .* Fp) .* Yr1);

    %% Vmat
    V_coef0_f = (1i .* nu_phi .* v ./ (root_u .* root_Delta .* denom_v)) * (mrow .* n2row) ...
                + (nu_v ./ (Delta.^2 .* denom_v)) .* (nrow .* (u.^2 .* (1 + nrow .* v.^2) - n2row .* v.^4)) ...
                - (nu_u .* root_u .* u .* v ./ (Delta.^2 .* root_v)) * nrow;

    V_coef0_fp = (nu_v .* u .* denom_u ./ Delta.^2) * ones_row ...
                + (nu_u .* root_u .* v ./ (Delta.^2 .* root_v)) .* (nrow .* Delta + denom_v);

    V_coef1_f = (1i .* nu_phi ./ (root_u .* root_Delta .* denom_v)) * (mrow .* alpha_row) ...
                - (nu_v .* v ./ (Delta.^2 .* denom_v)) .* alpha_row .* (-1 - 2 .* (nrow + 1) .* u.^2 + (2 .* nrow + 3) .* v.^2) ...
                - (nu_u .* root_u .* u ./ (Delta.^2 .* root_v)) * alpha_row;

    V_coef1_fp = (nu_u .* root_u ./ (Delta .* root_v)) * alpha_row;
    V_coef2_f = (nu_v ./ (Delta .* denom_v)) * (alpha_row .* beta_row);

    Vmat = common_coeffs .* ( ...
        (V_coef0_f .* Fr + V_coef0_fp .* Fp) .* Yr0 + ...
        (V_coef1_f .* Fr + V_coef1_fp .* Fp) .* Yr1 + ...
        (V_coef2_f .* Fr) .* Yr2);

    %% PHImat
    PHI_coef0_fp = (nu_phi .* u ./ Delta) * ones_row ...
                    + (1i .* nu_u ./ (root_v .* root_Delta)) * mrow;

    PHI_coef0_f = -(nu_phi ./ (denom_u .* denom_v)) * m2row ...
                - (nu_phi .* v.^2 ./ (denom_v .* Delta)) * nrow ...
                + (1i .* nu_v .* v ./ (root_u .* denom_v .* root_Delta)) * (mrow .* n2row) ...
                - (1i .* nu_u .* u ./ (denom_u .* root_v .* root_Delta)) * mrow;

    PHI_coef1_f = -(nu_phi .* v ./ (denom_v .* Delta)) * alpha_row ...
                + (1i .* nu_v ./ (root_u .* denom_v .* root_Delta)) * (mrow .* alpha_row);

    PHImat = common_coeffs .* ( ...
        (PHI_coef0_f .* Fr + PHI_coef0_fp .* Fp) .* Yr0 + ...
        (PHI_coef1_f .* Fr) .* Yr1);
end
end
