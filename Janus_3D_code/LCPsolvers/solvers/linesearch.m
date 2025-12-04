%% Nocedal and Wright pg 60
% Line search algorithm to satify the strong Wolfe-Conditions
function [x_k, Ax_k, f_k, grad_k, eta, p, Ap] = linesearch(x_km1, ...
    Ax_km1, f_km1, grad_km1, step, opts, debug)
if ~exist('debug', 'var') || isempty(debug)
    debug = false;
end
% When MVP's are less expensive a line search could be performed.
linesearchBudget = opts.linesearch.budget; c1 = opts.linesearch.c1; 
c2 = opts.linesearch.c2; tol = opts.linesearch.tol;
% phi the one dimensional restriction f(x_km1 + t*p(t))
% psi is the (relaxed) linearization f(x_km1) + c1*t*\nabla f(x_km1)^T*p(t)
b = opts.b;
function [phi, psi] = phi_and_psi(t)
    [x_k, Ax_k, f_k, grad_k, eta, p, Ap] = step(t, opts);
    phi = 1/2*dot(Ax_km1 + eta*Ap, x_km1 + eta*p) + dot(b, x_km1 + eta*p);
    psi =  f_km1 + c1*t*dot(grad_km1, p);
    % Update memory for eta prime estimator
    TT(end+1) = t; %#ok<*AGROW>
    PHI(end+1) = phi;
    PSI(end+1) = psi;
    [TT,ix] = sort(TT);
    PHI = PHI(ix);
    PSI = PSI(ix);
end
% The curvature condition is harder to check because eta and p are
% functions of t. So we approximate eta with a spline and p with
function ret = phiPrime(t)
    if numel(TT) >= 4
        phiFn = interp1(TT, PHI, 'pchip', 'pp');
    else
        phiFn = interp1(TT, PHI, 'linear', 'pp');
    end
    % Take the derivative of the spline
    phiFn.order=phiFn.order-1;
    phiFn.coefs=phiFn.coefs(:,1:end-1).*(phiFn.order:-1:1);
    ret = ppval(phiFn,t);
end
function ret = psiPrime(t)
    if numel(TT) >= 4
        psiFn = interp1(TT, PSI, 'pchip', 'pp');
    else
        psiFn = interp1(TT, PSI, 'linear', 'pp');
    end
    % Take the derivative of the spline
    psiFn.order=psiFn.order-1;
    psiFn.coefs=psiFn.coefs(:,1:end-1).*(psiFn.order:-1:1);
    ret = ppval(psiFn,t);
end
TT = 0;
PHI = f_km1;
PSI = f_km1;
tLo = 0;
phiLo = f_km1;
t = 1; % end of pg 59 of N&W
l = 1;

while l <= linesearchBudget
    [phi, psi] = phi_and_psi(t);
    if norm(p) < tol
        % The step is so small line search is pointless
        break
    elseif (phi > psi || ...
            (l > 1 && phiLo + tol < phi))
        zoom(t);
        break
    elseif abs(phiPrime(t)) <= c2/c1*abs(psiPrime(t))
        % Strong Wolf is satisfied
        % phi(t) < psi(t) && |phiPrime(t)| < |c2*phiPrime(0)|
        break
    elseif phiPrime(t) >= 0
        zoom(t);
        break
    end
    tLo = t;
    phiLo = phi;
    t = t*2;
    l = l+1;
end
if debug && l > linesearchBudget
    warning('Linesearch failed')
end % zoom

function zoom(tHi)
    l = l+1;
    while l <= linesearchBudget && (tHi - tLo) >= tol
        t = tLo + (tHi - tLo)/2;
        if l >= 4
            tt = linspace(tLo, tHi,100);
            phiInterp = interp1(TT, PHI, tt, 'pchip');
            [~, i] = min(phiInterp);
            if tt(i) > tLo+tol && tt(i) < tHi-tol
                t = tt(i);
            end
        end
        [phi, psi] = phi_and_psi(t);
        if norm(p) < tol
            % The step is so small line search is pointless
            return
        elseif (phi > psi || ...
                phi >= phiLo )
            % There is insufficient decrease or
            tHi = t;
        elseif abs(phiPrime(t)) <= c2/c1*abs(psiPrime(t))
            return
        elseif phiPrime(t) >= 0
            % In this regime the slope of phi is too positive.
            tHi = t;
        else
            % In this regime the slope of phi is too negative.
            tLo = t;
            phiLo = phi;
        end
        l = l+1;
    end
end % zoom

if debug && l >= linesearchBudget
    plotLineSearch(t, x_km1, Ax_km1, f_km1, grad_km1, step, opts)
end
end % lineSearch