%{
    Generates convergence plots and boxplots similar to the ones used in 
    an old paper.
%}
clear
DO_ON_SURFACE_PLOT = true;
DO_OFF_SURFACE_PLOT = false;
OBLATE_FLAG = true;

parr = [4, 6, 8, 10, 12, 14, 16];
distances = [3 2 10.^(0:-1:-3)];
% u0 = 100001/sqrt(200001); % AR = 1+1e-5
% u0 = 13/sqrt(69); % AR = 1.3
% u0 = 2/sqrt(3); % AR = 2
u0 = 4/sqrt(15); % AR = 4
aspect_ratio = round(u0 / sqrt(u0^2 - 1), 6);

%%%
%%% ON-SURFACE (comparing with kernelD)
%%%

if DO_ON_SURFACE_PLOT

% Find on-surface errors
errs_SPP = cell(1, numel(parr));
errs_surfdiv = cell(1, numel(parr));
for i = 1:numel(parr)
    p = parr(i);
    np = 2*p*(p+1);
    params = SpheroidalParameters;
    params.u0 = u0;
    params.a = 1/u0;
    params.oblate = OBLATE_FLAG; 
    params.centers = [0 0 0];
    params.isReal = true;

    if params.oblate
        X_self = oblate_spheroid_shape(p, params.u0, params.a);
    else
        X_self = prolate_spheroid_shape(p, params.u0, params.a);
    end
    nu_self = get_norm_vecs(p, params.u0, params.oblate);

    [u, v] = gl_grid(p);
    sigma_x = cos(u) .* sin(u).^2 + sin(u) .* cos(v);
    sigma_y = cos(u).^2 + sin(u).*cos(u).*sin(v);
    sigma_z = cos(u) .* sin(u).^2 - sin(u) .* cos(v) + exp(cos(u).^3);

    %%% Spectral calculation
    target_pts = cell(1, 1);
    target_pts{1} = X_self;
    target_nu = cell(1, 1);
    target_nu{1} = nu_self;
    [L2StkTLPx, L2StkTLPy, L2StkTLPz] = L2StkTLP(target_pts, target_nu, params, sigma_x, sigma_y, sigma_z, 1, false);
    [L2StkTLPdivx, L2StkTLPdivy, L2StkTLPdivz] = L2StkTLP(target_pts, target_nu, params, sigma_x, sigma_y, sigma_z, 1, true);

    %%% kerneldS calculation
    Df = kerneldS([], SurfaceSph(X_self));
    kerneldS_result = Df * reshape([sigma_x, sigma_y, sigma_z].', [], 1);
    kerneldS_result = reshape(kerneldS_result,3,[]).';
    
    errs_SPP{i} = norm(abs([L2StkTLPx{1}, L2StkTLPy{1}, L2StkTLPz{1}] - kerneldS_result))/norm(kerneldS_result);
    errs_surfdiv{i} = norm(abs([L2StkTLPdivx{1}, L2StkTLPdivy{1}, L2StkTLPdivz{1}] - kerneldS_result))/norm(kerneldS_result);
end

figure;
semilogy(parr, [errs_SPP{:}], 'o', parr, [errs_surfdiv{:}], 'o')
legend("Second-order derivative implementation", "Surface divergence implementation")
xlabel("p");
ylabel("relative error")
title("on-surface evaluation (w.r.t kerneldS) for aspect ratio="+string(aspect_ratio));
grid on;
hold off;

end % DO_ON_SURFACE_PLOT

%%%
%%% OFF-SURFACE (comparing with Kernel_Eval)
%%%

if DO_OFF_SURFACE_PLOT

parr = [16];
% Each cell contains a matrix corresponding to a p, where each column 
% represents the error associated with a specific distance from the
% surface.
errs_SPP_offsurf = cell(1, numel(parr));
errs_surfdiv_offsurf = cell(1, numel(parr));
for i = 1:numel(parr)
    p = parr(i);
    np = 2*p*(p+1);
    params = SpheroidalParameters;
    params.p = p;
    params.u0 = u0;
    params.a = 1/u0;
    params.oblate = OBLATE_FLAG; 
    params.centers = [0 0 0];

    [u, v] = gl_grid(p);
    sigma_x = cos(u) .* sin(u).^2 + sin(u) .* cos(v);
    sigma_y = cos(u).^2 + sin(u).*cos(u).*sin(v);
    sigma_z = cos(u) .* sin(u).^2 - sin(u) .* cos(v) + exp(cos(u).^3);

    nu_trg = get_norm_vecs(p, params.u0, params.oblate);
    errs_SPP_offsurf{i} = cell(1, numel(distances));
    errs_surfdiv_offsurf{i} = cell(1, numel(distances));
    for j = 1:numel(distances)
        d = distances(j);
        if params.oblate
            X_trg = oblate_spheroid_shape(p, params.u0, params.a);
        else
            X_trg = prolate_spheroid_shape(p, params.u0, params.a);
        end
        X_trg = X_trg + d*nu_trg;
        nt = size(X_trg, 1);
    
        %%% Stokes Kernel_Eval
        % Get source geometry and weights
        if params.oblate
            X_src_orig = oblate_spheroid_shape(p, params.u0, params.a);
        else
            X_src_orig = prolate_spheroid_shape(p, params.u0, params.a);
        end
        N_src_orig = params.get_Norm(p, 1);
    
        % Quadrature weights for Kernel_Eval
        Sns = SurfaceSph(X_src_orig);
        [~, gwt_gl] = g_grid(p + 1);
        wt_gl = pi/p * repmat(gwt_gl', 2*p, 1) ./ sin(gl_grid(p));
        wt_gl = wt_gl(:);
        W_src_orig = Sns.geoProp.W .* wt_gl;
        Wv = repmat(W_src_orig,1,3)'; Wv=Wv(:);
    
        pot = 'TSL_Stk_3D';
        KEparams = Kernel_Eval_parameters(pot,0,1,1,1,1e-12,2,400,1);
        KEparams.dim = 3;
        Xv = reshape(repmat(X_src_orig,1,3)',3,[])';
        KEparams.X = Xv;
        KEparams.W2 = Wv.';
        KEparams.nor = reshape(repmat(nu_trg,1,3)',3,[])';
        KEparams.ci = repmat((1:3)', size(X_trg, 1), 1);
        KEparams.cj = repmat((1:3)', size(X_src_orig, 1), 1);
    
        Xtrg_ii = reshape(repmat(X_trg,1,3)',3,[])';
        TLP_mat = Kernel_Eval(Xtrg_ii, Xv, KEparams);
    
        %%% Compute results
        sig = reshape([sigma_x,sigma_y,sigma_z].',[],1);
        TLP_kernel_eval = TLP_mat * sig;
        TLP_kernel_eval = reshape(TLP_kernel_eval,3,[]).';
    
        %%% Spectral calculation
        target_pts = cell(1, 1);
        target_pts{1} = X_trg;
        target_nu = cell(1, 1);
        target_nu{1} = nu_trg;
        [L2StkTLPx, L2StkTLPy, L2StkTLPz] = L2StkTLP(target_pts, target_nu, params, sigma_x, sigma_y, sigma_z, 1, false);
        [L2StkTLPdivx, L2StkTLPdivy, L2StkTLPdivz] = L2StkTLP(target_pts, target_nu, params, sigma_x, sigma_y, sigma_z, 1, true);
        
        %%% Compare results
        errs_SPP_offsurf{i}{j} = log10(abs([L2StkTLPx{1}, L2StkTLPy{1}, L2StkTLPz{1}] - TLP_kernel_eval) ./ abs(TLP_kernel_eval));
        errs_surfdiv_offsurf{i}{j} = log10(abs([L2StkTLPdivx{1}, L2StkTLPdivy{1}, L2StkTLPdivz{1}] - TLP_kernel_eval) ./ abs(TLP_kernel_eval));
    end
end

%%% Plot
for i = 1:numel(parr)
    p_val = parr(i);
    
    % Pre-calculate Y-axis limits
    min_y = inf;
    max_y = -inf;
    for j = 1:numel(distances)
        all_data_at_dist = [errs_SPP_offsurf{i}{j}(:); errs_surfdiv_offsurf{i}{j}(:)];
        if ~isempty(all_data_at_dist)
            min_y = min(min_y, min(all_data_at_dist));
            max_y = max(max_y, max(all_data_at_dist));
        end
    end
    padding = (max_y - min_y) * 0.05;
    y_limits = [min_y - padding, max_y + padding];

    % Create a figure for the current p-value
    figure('Name', sprintf('p = %d', p_val), 'Position', [100, 100, 1200, 500]);
    sgtitle(sprintf('off-surface evaluation with p=%d, aspect ratio='+string(aspect_ratio), p_val));
    
    num_plots = numel(distances);
    
    % Loop through distances to create a subplot for each
    for j = 1:numel(distances)
        d = distances(j);
        subplot(1, num_plots, j);
        
        data1 = errs_SPP_offsurf{i}{j}(:);
        data2 = errs_surfdiv_offsurf{i}{j}(:);
        combined_data = [data1; data2];
        grouping_variable = [
            repmat({'Second-order derivative'}, numel(data1), 1) ;
            repmat({'Surface divergence implementation'}, numel(data2), 1)
        ];
        
        % Create the box chart
        boxchart(ones(size(combined_data)), combined_data, 'GroupByColor', categorical(grouping_variable));
        
        % --- Formatting ---
        xticklabels('')
        ylim(y_limits);
        title(sprintf('d=%.0d', d));
        grid on;

        if j == 1
            ylabel('log10(relative error)');
        end
    end

    legend
end

end % DO_OFF_SURFACE_PLOT