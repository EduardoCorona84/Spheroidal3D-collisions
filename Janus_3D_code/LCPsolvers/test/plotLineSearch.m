function plotLineSearch(tstar, x_km1, Ax_km1, f_km1, grad_km1, step, opts)
c1 = opts.linesearch.c1; 
c2 = opts.linesearch.c2; 
% phi the one dimensional restriction f(x_km1 + t*p(t))
% psi is the (relaxed) linearization f(x_km1) + c1*t*\nabla f(x_km1)^T*p(t)
function [phi, phiHat, psi] = phi_and_psi(t)
    b = opts.b;
    opts.stepSize.bwd = 'uniform';
    [x_k, Ax_k, f_k, grad_k, eta, p, Ap] = step(t, opts);
    phi = 1/2*dot(Ax_km1 + eta*Ap, x_km1 + eta*p) + dot(b, x_km1 + eta*p);
    opts.stepSize.bwd = 'opt';
    [x_k, Ax_k, f_k, grad_k, eta, p, Ap] = step(t, opts);
    phiHat = 1/2*dot(Ax_km1 + eta*Ap, x_km1 + eta*p) + dot(b, x_km1 + eta*p);
    psi =  f_km1 + c1*t*dot(grad_km1, p);
    % Update memory for eta prime estimator
    TT(end+1) = t; %#ok<*AGROW>
    PHI(end+1) = phi;
    PHI_HAT(end+1) = phiHat;
    PSI(end+1) = psi;
    [TT,ix] = sort(TT);
    PHI = PHI(ix);
    PHI_HAT = PHI_HAT(ix);
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

function ret = phiHatPrime(t)
    if numel(TT) >= 4
        phiFn = interp1(TT, PHI_HAT, 'pchip', 'pp');
    else
        phiFn = interp1(TT, PHI_HAT, 'linear', 'pp');
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
PHI_HAT = f_km1;
PSI = f_km1;
%% Plot Curvilinear linesearch 
tt = [10.0.^(-8:1:-3) 0.01:0.01:10, linspace(tstar-tstar,tstar+tstar, 100) tstar];
tt = sort(unique(tt));
tt = tt(tt>0);
goodtt = [];
goodPhi =[];
goodttHat = [];
goodPhiHat = [];
phiBest = Inf;
for i = 1:numel(tt)
    t = tt(i);
    %% Don't look for a beter step size after prox
    [phi, phiHat, psi] = phi_and_psi(t);
    if abs(phiPrime(t)) <= c2/c1*abs(psiPrime(t))
        goodtt(end+1) = t;
        goodPhi(end+1) = phi;
    end
    if abs(phiHatPrime(t)) <= c2/c1*abs(psiPrime(t))
        goodttHat(end+1) = t;
        goodPhiHat(end+1) = phiHat;
    end
    if phi < phiBest 
        phiBest = phi;
        tBest = t; 
    end
    if t == tstar
        phiUs = phiHat;
    end
    if phi > psi && t > tstar
        break
    end
end
%% Plot
f = gcf;
set(f,"Position",[0,0,800,900])
clf;
subF =subplot('Position',[0.15, 0.12,.9,.78]);
hold on
ps = {};
% ps{end+1} = plot(tt, psiq, 'LineWidth',10,'Color','#ff7f0e');
ps{end+1} = plot(TT, PSI, 'LineWidth',10,'Color','#d62728');
ps{end+1} = plot(TT, PHI, 'LineWidth',10, 'Color','#1f77b4');
ps{end+1} = plot(TT, PHI_HAT,'LineWidth',10, 'Color','#17becf');
ps{end+1} = plot(goodtt, goodPhi, 'Color', '#565A5C', 'LineWidth',10, 'LineStyle', ':');
ps{end+1} = plot(goodttHat, goodPhiHat, 'Color', '#565A5C', 'LineWidth',10, 'LineStyle', ':');
% Markers
ps{end+1} = plot(tBest, phiBest, 'LineStyle','none', 'Marker','square', 'MarkerSize',40,...
    'Color','black', 'MarkerFaceColor','black');
ps{end+1} = plot(tstar, phiUs, 'LineStyle','none', 'Marker','diamond', 'MarkerSize',40, ...
    'Color','#2ca02c', 'MarkerFaceColor','#2ca02c');
a = min([PSI PHI PHI_HAT]);
b = max([PSI PHI PHI_HAT]);
d = b-a;
% a = max([1e-2, a-.1*p]);
a = a-.1*d;
b = b +.1*d;
ylim([a,b])
ax = gca;
ax.FontSize = 35;
xlabel('$t$', 'FontSize', 50, 'Interpreter','latex');
ylabel('', 'FontSize', 50, 'Interpreter','latex');
ytickformat('%.3g')
xlim([TT(1), TT(end)])
set(subF, 'Position', [0.15, 0.12,.8,.725])
title(['$\phi(t^*)$ = ' sprintf('%.5g', phiBest) ...
    ', $\hat{\phi}(\hat{t})$ = ' sprintf('%.5g', phiUs)], 'FontSize', 50, 'Interpreter','latex');
k = opts.k;
sgtitle(sprintf("Iteration %d", k), 'FontSize',50);

% saveas(f, sprintf('/Users/niru8088/scratch/Spheroidal3D-collisions/docs/fig/linesearch_%p.png', k))

% if k ==1
%     for i = 1:numel(ps)
%         ps{i}.XData = NaN;
%         ps{i}.YData = NaN;
%     end
%     set(gca, 'visible', 'off')
%     sgtitle('')
%     legend( ...
%         { ...
%         % '$\psi_{\mathbf{q}^{(k)}}(t)$'; ...
%         '$\psi(t)$'; ...
%         '$\phi(t)$'; ...
%         '$\hat{\phi}(t)$';  ...
%         '$|\phi^\prime(t)| < \frac{c_2}{c_1} |\psi^\prime(t)|$'; ...
%         '$|\hat{\phi}^\prime(t)| < \frac{c_2}{c_1} |\psi^\prime(t)|$'; ...
%         '$\phi(t^*)$'; ...
%         % '$\hat{\phi}(t^*)$'; ...
%         '$\hat{\phi}(1)$' ...
%         }, ...
%         'FontSize', 50, 'Interpreter','latex','Position',[0.1, 0.1,.8,.8], ...
%         'TextColor','black', 'Color','white')
%     % saveas(f, sprintf('/Users/niru8088/scratch/Spheroidal3D-collisions/docs/fig/linesearch_legend.png'))
% end
end % plotLineSearch
