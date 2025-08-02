function [kappa, eval, rel_residual_step, abs_residual_step, info] = stepSize(k, p, grad_k, opts, info, iter, s_k, y_k, tol_type)

if k == 0
    mode = opts.kappa.init;
elseif k < 0
    mode = opts.kappa.bwd;
else
    mode = opts.kappa.fwd;
end

if k > 0 && (~exist('s_k','var') || ~exist('y_k','var'))
    assert(~contains(lower(mode),'bb'), ['BB steps require s_k and y_k '...
        'which are not available when k == 0']);
end

switch lower(mode)
    case 'bb1'
        kappa = (s_k'*s_k)/(s_k'*y_k);
    case 'bb2'
        kappa = (s_k'*y_k)/(y_k'*y_k);
    case 'opt'
        % For QP, this is the optimal step length (see page 56 of n&W)
        [eval, gmres_iters_iter, rel_residual_step, abs_residual_step] = opts.A(p, opts.noise_control.noise, tol_type);
        info.gmres_iters(iter + 1) = info.gmres_iters(iter + 1) + gmres_iters_iter;
        if dot(p, eval) < 0
            kappa = (s_k'*s_k)/(s_k'*y_k);
        else
            kappa = -dot(p, grad_k) / dot(p, eval);
        end
    case 'uniform' 
        kappa = 1;
    otherwise 
        % NIC: Taken from projected gradient descent... I think this a bug?
        % [~,ytmp] = fg(y);
        % kappa = (y'*y) / (y'*ytmp);
        % NIC: In general for the first iteration without Lipschitz info, we just
        % use 1?
        error([mode ' is not a valid step size rule'])
end