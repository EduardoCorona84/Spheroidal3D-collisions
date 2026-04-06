function [Akk, Tself0] = SpheroidalMS_BuildTDSelfBlock(parbd, Lk, body_idx, Tself0)
%{
Build the current-frame TD self block for one spheroidal body.

The TSL self-interaction is assembled once in the body's local
frame and can be passed back in through Tself0 for reuse. The current
self block is then obtained by rotating that local block into the
body's current frame and adding 0.5I + Lk.

The key idea is that the TSL operator, as used by spheroidal_mobility, is
in the global frame, so we need to account for this.
%}

Nb = parbd.Nb;
np = parbd.np;

if nargin < 4 || isempty(Tself0)
    Tself0 = LOCAL_build_td_local_tself(parbd, body_idx);
end

Qk = parbd.MRot{body_idx};
if isequal(Qk, eye(3))
    Tself = Tself0;
else
    Tself = Rotate_Operator(Tself0, Qk, np);
end

idx = (1:Nb) + Nb * (body_idx - 1);
Akk = Tself + 0.5 * eye(Nb) + Lk(idx, idx);
end %% END MAIN FUNCTION

function Tself = LOCAL_build_td_local_tself(parbd, body_idx)
    np = parbd.np;
    Nb = parbd.Nb;

    if isfield(parbd, 'tsl_dealiasing')
        tsl_dealiasing_flag = parbd.tsl_dealiasing;
    else
        tsl_dealiasing_flag = true;
    end
    if isfield(parbd, 'tsl_dealiasing_pad')
        tsl_dealiasing_pad = parbd.tsl_dealiasing_pad;
    else
        tsl_dealiasing_pad = 4;
    end
    if isfield(parbd, 'tsl_backend') && ~isempty(parbd.tsl_backend)
        tsl_backend = parbd.tsl_backend;
    else
        tsl_backend = 'cartesian';
    end

    prm = zeros(1, Nb);
    prm(1:np) = 1:3:Nb;
    prm(np+1:2*np) = 2:3:Nb;
    prm(2*np+1:3*np) = 3:3:Nb;
    iprm = zeros(1, Nb);
    iprm(prm) = 1:Nb;

    shape_type = parbd.shape_type(body_idx);
    if strcmp(shape_type, 'sphere')
        error('Self-block build currently supports prolate/oblate only.');
    end

    [u0, a] = calculate_u0_and_a_from_radii(shape_type, parbd.equ_radii(body_idx), parbd.polar_radii(body_idx));
    oblate = strcmp(shape_type, 'oblate');

    params_i = SpheroidalParameters();
    params_i.sigma = eye(np);
    params_i.u0 = u0;
    params_i.a = a;
    params_i.oblate = oblate;
    params_i.centers = [0 0 0];
    params_i.thetas = 0;
    params_i.phis = 0;
    params_i.Rmat = eye(3);

    I = eye(np);
    Z = zeros(np);
    source_Gmatrix = sparse(Gmatrix(params_i.p, u0, 0, oblate));

    Yx = LOCAL_eval_l2stk_dense_self(params_i, I, Z, Z, source_Gmatrix, np, tsl_backend, tsl_dealiasing_flag, tsl_dealiasing_pad);
    Yy = LOCAL_eval_l2stk_dense_self(params_i, Z, I, Z, source_Gmatrix, np, tsl_backend, tsl_dealiasing_flag, tsl_dealiasing_pad);
    Yz = LOCAL_eval_l2stk_dense_self(params_i, Z, Z, I, source_Gmatrix, np, tsl_backend, tsl_dealiasing_flag, tsl_dealiasing_pad);

    Tself = [Yx, Yy, Yz];
    Tself = Tself(:, iprm);
end

function Yblk = LOCAL_eval_l2stk_dense_self(params_i, sigma_x, sigma_y, sigma_z, source_Gmatrix, np, tsl_backend, tsl_dealiasing_flag, tsl_dealiasing_pad)
    switch tsl_backend
        case 'spheroidal'
            [vx, vy, vz] = ...
                L2StkTLPOptimized([], [], params_i, sigma_x, sigma_y, sigma_z, source_Gmatrix, false, tsl_dealiasing_flag, tsl_dealiasing_pad);
        case 'cartesian'
            [vx, vy, vz] = ...
                L2StkTLPCartesianOptimized([], [], params_i, sigma_x, sigma_y, sigma_z, source_Gmatrix, tsl_dealiasing_flag, tsl_dealiasing_pad);
        otherwise
            error('Unknown TSL backend "%s".', tsl_backend);
    end
    y = reshape([vx(:), vy(:), vz(:)].', [], 1);
    Yblk = reshape(y, 3*np, np);
end
