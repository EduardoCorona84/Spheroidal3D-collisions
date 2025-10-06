function [kappa, Ap] = stepSize(k, p, x, Ax, opts, s_k, y_k, Ap)

if k == 1
    mode = opts.stepSize.init;
elseif k > 1 %  fwd
    mode = opts.stepSize.kappa;
else % bwd
    mode = opts.stepSize.eta;
end

if k > 0 && (~exist('s_k','var') || ~exist('y_k','var'))
    assert(~contains(lower(mode),'bb'), ['BB steps require s_k and y_k '...
        'which are not available when k == 0']);
end
if ~exist('Ap','var')
    Ap = [];
end
switch lower(mode)
    case 'bb1'
        kappa = (s_k'*s_k)/(s_k'*y_k);
    case 'bb2'
        kappa = (s_k'*y_k)/(y_k'*y_k);
    case 'opt'
        % For QP, this is the optsimal step length (see page 56 of n&W)
        % Notice that is A is not perfectly symmetric, then we do not have
        % kappa = -(Ax+b)'p/(p'Ap)
        Ap = opts.A(p);
        kappa = -(1/2*(dot(p, Ax) +  dot(x, Ap)) + dot(p,opts.b)) / dot(p, Ap);
        % In the bwd case, we need to stay in the feasible set
        % - Because x>0 and x + p > 0, via convexity kappa \in (0,1] is good
        % - In the other case, we need to check when the ray intersects the
        %   positive orthant this is separable, and we can find when each element
        %   of x + kappa p = 0 by taking -x / p elementwize. When -x / p < 0 then
        %   it is irrelevant. But if not then we need to make sure that we only
        %   travel to the closest feasible point.
        if k < 0
            if kappa <= 1
                return
            else
                kappa_list = - x(p < 0) ./ p(p < 0);
                if isempty(kappa_list)
                    return;
                end
                kappa = min(kappa, min(kappa_list));
            end
            
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
end