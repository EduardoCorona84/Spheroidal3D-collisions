function [t, Ap] = stepSize(k, p, x, Ax, opts, s, y)

if k > 0 %  fwd
    mode = opts.stepSize.fwd;
else % bwd
    mode = opts.stepSize.bwd;
end

if k > 0 && (~exist('s','var') || ~exist('y','var'))
    assert(~contains(lower(mode),'bb'), ['BB steps require s and y '...
        'which are not available when k == 0']);
end

Ap = [];
switch lower(mode)
    case 'bb1'
        if isempty(s)
            t = 1;
            return 
        end 
        t = (s'*s)/(s'*y);
    case 'bb2'
        if isempty(s)
            t = 1;
            return 
        end
        t = (s'*y)/(y'*y);
    case 'lipschitz'
        t = 1/opts.stepSize.L;
    case 'opt'
        % For QP, this is the optsimal step length (see page 56 of n&W)
        % Notice that is A is not perfectly symmetric, then we do not have
        % t = -(Ax+b)'p/(p'Ap)
        b = opts.b;
        if dot(p, Ax+b) > -1e-8
            % p neads to be a descent direction
            % if it is not sufficiently collinear with the gradient then
            % return 1
            t = 1; 
        else
            Ap = opts.A(p);
            t = -(dot(p, Ax+b)) / dot(p, Ap);
            % In the bwd case, we need to stay in the feasible set
            % - Because x>0 and x + p > 0, via convexity t \in (0,1] is good
            % - In the other case, we need to check when the ray intersects the
            %   positive orthant this is separable, and we can find when each element
            %   of x + t p = 0 by taking -x / p elementwize. When -x / p < 0 then
            %   it is irrelevant. But if not then we need to make sure that we only
            %   travel to the closest feasible point.
            if k < 0
                mask = p < 0;
                if ~any(mask)
                    return;
                end
                t = min([t; - x(mask) ./ p(mask)]);
            end
        end
    case 'uniform' 
        t = 1;
    otherwise 
        error([mode ' is not a valid step size rule'])
end
if nargout == 2 && isempty(Ap)
    Ap = opts.A(p);
end

end