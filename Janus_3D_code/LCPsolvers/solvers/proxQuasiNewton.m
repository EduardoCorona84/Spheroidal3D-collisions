function [x, info, opts] = proxQuasiNewton(fg, x0, opts)
[opts, info] = defaultLCPOpts(opts, x0);
checkOpts(opts)
n = numel(x0); k = 0;
x_k = x0; Ax_k = zeros(n,1); s = []; y = [];
if ~all(x0 == 0) && ~isfield(opts, 'Ax_k') % No iterate history
    Ax_k = opts.A(x_k);
    s = x_k;
    y = Ax_k;
elseif isfield(opts, 'Ax_k') 
    % Has iterate history, assuming solving a subproblem
    % No secant condition to update as it is stored from the previous subproblem.
    Ax_k = opts.Ax_k;
end
[f_k, grad_k] = fg(x_k, Ax_k);
while true
    [converged, info] = checkConvergence(k, f_k, x_k, ...
        grad_k, [], info, opts);
    if converged
        x = x_k; 
        % Store the last iterate for the next subproblem and update the memory. 
        opts.Ax_k = Ax_k;
        [~, ~, ~, ~, opts] = updateQNMemory(s, y, opts);
        break
    end
    % Increase k 
    k = k + 1;
    opts.k = k;
    x_km1 = x_k; Ax_km1 = Ax_k; f_km1= f_k; grad_km1 = grad_k;
    % Update (inverse) Hessian Approximation
    [H, opts] = updateHk(s, y, opts);
    % Define this step 
    q = -H(grad_km1);
    prox_k = @(xtilde) prox(xtilde, opts);
    step_k = @(t, opts) fwdBwdstep(t, x_km1, Ax_km1, q, prox_k, fg, opts);
    % Linesearch 
    [x_k, Ax_k, f_k, grad_k] = linesearch(x_km1, Ax_km1, f_km1, grad_km1, ...
        step_k, opts);
    % Save secant conditions
    s = x_k - x_km1;
    y = grad_k - grad_km1;
end
end % proxQuasiNewton

function checkOpts(opts) %#ok<INUSD>

end