function [res_x, res_y, res_z] = L2StkSLPOptimized(X_trg, params, sigma_x, sigma_y, sigma_z, source_Gmatrix)
%{
Optimized version of L2Stk for the Stokes single layer potential.
%}

p = params.p;
u0 = params.u0;
a = params.a;
isReal = params.isReal;
oblate = params.oblate;

[np, nf, ns] = size(sigma_x);

if ~isempty(X_trg) && (size(X_trg, 2) ~= 3 || size(X_trg, 3) ~= ns)
    error("X_trg must be nt x 3 x ns.");
end

% Compute y dot sigma on the source surface (local frame)
y_dot_sig = zeros(np, nf, ns);
for i = 1:ns
    if params.oblate(i)
        Xloc_i = oblate_spheroid_shape(p, params.u0(i), params.a(i));
    else
        Xloc_i = prolate_spheroid_shape(p, params.u0(i), params.a(i));
    end
    y_dot_sig(:,:,i) = sigma_x(:,:,i).*Xloc_i(:,1) + sigma_y(:,:,i).*Xloc_i(:,2) + sigma_z(:,:,i).*Xloc_i(:,3);
end

% Compute SH coefficients and Gmatrix transforms
[shc_x, shc_y, shc_z, shc_ydotsig] = generate_shc(sigma_x, sigma_y, sigma_z, y_dot_sig);
[Gshc_x, Gshc_y, Gshc_z, Gshc_ydotsig] = generate_mod_shc(source_Gmatrix, shc_x, shc_y, shc_z, shc_ydotsig);

% Combine the densities so we can evaluate in one pass
Gshc_all = cat(2, Gshc_x, Gshc_y, Gshc_z, Gshc_ydotsig);

% Build unit vectors for SP
[nu_x, nu_y, nu_z] = build_unit_nu(X_trg, np, ns);

% Evaluate SL and SP in one pass each
SL_all = spheroidalSLOptimized(p, u0, a, oblate, Gshc_all, isReal, X_trg);
[SPx_all, SPy_all, SPz_all] = spheroidalSPOptimized(p, u0, a, oblate, Gshc_all, isReal, nu_x, nu_y, nu_z, X_trg);

% Retrieve back the layer potentials
idx1 = 1:nf;
idx2 = nf + (1:nf);
idx3 = 2*nf + (1:nf);
idx4 = 3*nf + (1:nf);

SL1 = slice_density(SL_all, idx1);
SL2 = slice_density(SL_all, idx2);
SL3 = slice_density(SL_all, idx3);

SP1dx = slice_density(SPx_all, idx1);
SP2dx = slice_density(SPx_all, idx2);
SP3dx = slice_density(SPx_all, idx3);
Fdx   = slice_density(SPx_all, idx4);

SP1dy = slice_density(SPy_all, idx1);
SP2dy = slice_density(SPy_all, idx2);
SP3dy = slice_density(SPy_all, idx3);
Fdy   = slice_density(SPy_all, idx4);

SP1dz = slice_density(SPz_all, idx1);
SP2dz = slice_density(SPz_all, idx2);
SP3dz = slice_density(SPz_all, idx3);
Fdz   = slice_density(SPz_all, idx4);

% Assemble result
if isempty(X_trg)
    [Xloc, ~] = params.get_X();
    np_div = size(Xloc, 1) / ns;
    np_check = round(np_div);
    if abs(np_div - np_check) > 1e-8
        error("\n size of target on surface is not multiple of ns\n")
    end

    res_x = zeros(np_check, nf, ns);
    res_y = zeros(np_check, nf, ns);
    res_z = zeros(np_check, nf, ns);

    for i = 1:ns
        Xi = Xloc((i-1)*np_check+1:i*np_check, :);
        res_x(:,:,i) = 0.5 .* (SL1(:,:,i) - Xi(:,1).*SP1dx(:,:,i) - Xi(:,2).*SP2dx(:,:,i) - Xi(:,3).*SP3dx(:,:,i) + Fdx(:,:,i));
        res_y(:,:,i) = 0.5 .* (SL2(:,:,i) - Xi(:,1).*SP1dy(:,:,i) - Xi(:,2).*SP2dy(:,:,i) - Xi(:,3).*SP3dy(:,:,i) + Fdy(:,:,i));
        res_z(:,:,i) = 0.5 .* (SL3(:,:,i) - Xi(:,1).*SP1dz(:,:,i) - Xi(:,2).*SP2dz(:,:,i) - Xi(:,3).*SP3dz(:,:,i) + Fdz(:,:,i));
    end
else
    nt = size(X_trg, 1);
    res_x = zeros(nt, nf, ns);
    res_y = zeros(nt, nf, ns);
    res_z = zeros(nt, nf, ns);

    for i = 1:ns
        Xi = X_trg(:,:,i);
        res_x(:,:,i) = 0.5 .* (SL1(:,:,i) - Xi(:,1).*SP1dx(:,:,i) - Xi(:,2).*SP2dx(:,:,i) - Xi(:,3).*SP3dx(:,:,i) + Fdx(:,:,i));
        res_y(:,:,i) = 0.5 .* (SL2(:,:,i) - Xi(:,1).*SP1dy(:,:,i) - Xi(:,2).*SP2dy(:,:,i) - Xi(:,3).*SP3dy(:,:,i) + Fdy(:,:,i));
        res_z(:,:,i) = 0.5 .* (SL3(:,:,i) - Xi(:,1).*SP1dz(:,:,i) - Xi(:,2).*SP2dz(:,:,i) - Xi(:,3).*SP3dz(:,:,i) + Fdz(:,:,i));
    end
end
end %% END MAIN FUNCTION