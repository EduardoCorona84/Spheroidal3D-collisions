function xstar = prox_rank1(y, h0, w, sigma, opts, plotFlag)

if all(y >= 0)
    xstar = y;
    return 
elseif isempty(w)
     xstar = max(y, 0);
     return 
end

if ~exist('opts','var') || isempty(opts) || ~isfield(opts, 'prox')
    opts = struct('prox', struct('debug', false));
end

if ~exist('plotFlag','var') || isempty(plotFlag) 
    plotFlag = false;
end
%% Input checks
assert(size(y,2) == 1, 'y should be a column vector');
assert(size(h0,1) == 1 && size(h0,2) ==1, 'h0 should be a scalar');
assert(size(w,2) == 1, 'w should be a column vector');
assert(sigma == 1 || sigma == -1, 'Sigma should be the sign of the rank 1 update')
%% Find points of discontinuituy in the piecewise linear function
alphas = - sigma * y ./ (h0 .* w);
alphas = sort(alphas); % perhaps one day do something smarter than sorting...
%% Specify the values of the piecewise linear function $L$ at each value of
% alpha
la = @(a) y + sigma * a * h0 .* w;
N = length(y);
Ls = zeros(N,1);
for n = 1:N
    a = alphas(n);
    l = la(a);
    mask = l > 0;
    Ls(n) = a + dot(w, y) - dot(w(mask), l(mask));
end
%% Specify the interval where the root occurs
[ix, ~] = find(Ls < 0, 1, 'last');
if isempty(ix) % This means that all the constraints are active.
    a = alphas(1) - 10; % arbitrary
else
    a = alphas(ix);
end
l = la(a);
mask = l > 10*eps() | (abs(l) < 10*eps() & sigma * h0 * w > 0);
alphastar = (dot(w,y) - dot(w(mask), y(mask))) ...
    / (sigma * dot(w(mask), h0 * w(mask)) - 1 );

xstar = max(y + sigma * alphastar*h0.*w, 0);

if opts.prox.runCVX 
    %% idiot check
    % H = h0 + sigma * u * u'
    % B = 1 / h0 * I - sigma * w * w';
    n = length(y);
    B = eye(n) ./ h0 - sigma * w * w';
    try %#ok<TRYNC>
        cvx_begin quiet
            variable z(N)
            minimize( dot(y-z,B*(y-z)) )
            subject to 
                z >= 0
        cvx_end
        cvxErr = norm(z - xstar) / norm(xstar);
        assert(cvxErr < 1e-4, 'Our solution does not match CVX')
    end 
end 

if plotFlag
    amin = min([alphas; alphastar]); 
    amin = amin - abs(amin)*0.1;
    amax = max([alphas; alphastar]);
    amax = amax + abs(amax)*0.1;
    aa = sort([linspace(amin, amax, 1000) alphas']);
    LL = zeros(numel(aa),1);
    for n = 1:numel(aa)
        a = aa(n);
        l = la(a);
        mask = l > 0;
        LL(n) = a + dot(w, y) - dot(w(mask), l(mask));
    end
    plotlyjs_colors = {"#1f77b4";  % muted blue
        "#ff7f0e";  % safety orange
        "#2ca02c";  % cooked asparagus green
        "#d62728";  % brick red
        "#9467bd";  % muted purple
        "#8c564b";  % chestnut brown
        "#e377c2";  % raspberry yogurt pink
        "#7f7f7f";  % middle gray
        "#bcbd22";  % curry yellow-green
        "#17becf"};   % blue-teal
    gcf
    clf
    hold on
    plot(aa, zeros(numel(aa),1), '--', 'Color', '#565A5C', 'LineWidth',4)
    plot(aa, LL, 'Color', plotlyjs_colors{1}, 'LineWidth',6)
    plot(alphas, Ls, '.', 'Color', plotlyjs_colors{2}, 'MarkerSize', 40)
    plot(alphastar, 0, 'o', 'Color', plotlyjs_colors{3}, 'MarkerFaceColor', plotlyjs_colors{3}, 'MarkerSize', 15)
    legend({'0', '$\mathcal{L}(\alpha)$', '$\{\alpha_i\}$', '$\alpha^*$'},...
        'Interpreter','latex','Location','northwest','FontSize',30);

end

end % prox_rank1
