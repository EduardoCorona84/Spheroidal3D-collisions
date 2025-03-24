function y = prox_rank1_nn(x,d,u, plotFlag)
assert(all(d>=0), 'this code assumes that all di > 0');
if all(x > 0)
    y = x; 
    return 
end

if ~exist('plotFlag', 'var') || isempty(plotFlag)
    plotFlag = false; 
end

N = length(x); 
alphas = d .* x ./ u; 
alphasSorted = sort(alphas);
f = ones(N,1);
for i = 1:N
    ai = alphasSorted(i);
    % get all the alphas to the left of the current, grab all that are
    % decreasing -> xi - alpha*ui/di < 0
    % similar logic for those to the right
    tildeMask = (alphas <= ai & u >= 0) | ...
        (alphas >= ai & u <= 0);
    xtilde = x(tildeMask);
    utilde = u(tildeMask);

    hatMask = ~tildeMask;
    dhat = d(hatMask);
    uhat = u(hatMask);
    f(i) = dot(utilde, xtilde) + ...
        (1 + dot(uhat, uhat ./ dhat)) * ai;
    if f(i) > 0 
        break 
    end
end

j = find(f < 0, 1,'last');
aj = alphasSorted(1) - 1;
if ~isempty(j)
    aj = alphasSorted(j);
end
tildeMask = (alphas <= aj & u >= 0) | ...
    (alphas >= aj & u <= 0);
xtilde = x(tildeMask);
utilde = u(tildeMask);

hatMask = ~tildeMask;
dhat = d(hatMask);
uhat = u(hatMask);
alphastar =  -dot(utilde, xtilde) / (1 + dot(uhat, uhat ./ dhat));
fstar = (1 + dot(uhat, uhat ./ dhat)) * alphastar + dot(utilde, xtilde);

y = max(x - alphastar * u ./ d,0) ;

if plotFlag
    amin = min([alphas; alphastar]); 
    amin = amin - abs(amin)*0.1;
    amax = max([alphas; alphastar]);
    amax = amax + abs(amax)*0.1;
    aa = amin:.1:amax;
    fHandle = @(a) a + dot(u, x - max(0, x - a*u./d));
    ff = arrayfun(fHandle, aa);
    gcf
    clf
    hold on
    plot(aa,ff)
    plot(alphasSorted, f, '*')
    plot(alphastar, fstar, '*')
    legend({'$\mathcal{L}(\alpha)$', '$\{\alpha_i\}$', '$\alpha^*$'},...
        'Interpreter','latex');
end