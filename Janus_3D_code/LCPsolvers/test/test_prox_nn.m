plotFlag = true;
Ns = arrayfun(@(e) round(10^e), 1:0.5:3);
times = zeros(length(Ns));
for i = 1:length(Ns)
    % generate test problem
    N = Ns(i);
    x = 100*randn(N,1);
    d = abs(100*rand()) * ones(N,1);
    u = 100*randn(N,1);
    %Nic
    tic
    ynic = prox(x, d, u);
    times(i) = toc();
    % Stephen
    tic
    ystephen = proj_rank1_box( zeros(N,1), Inf*ones(N,1), x, d, u);
    stephentime = toc();
    % CVX
    b = 1 ./ d;
    denom = sqrt(1 + dot(u .* b, u));
    v = b .* u / denom;
    L = @(z) 1/2*dot(z-x, (diag(b) - v*v')*(z-x));
    tic
    cvx_begin quiet
        variable yCvx(N)
        minimize L(yCvx)
        subject to 
            0 <= yCvx

    cvx_end
    t2 = toc;
    % Display
    disp(['N = ' num2str(N)])
    disp(['  nic     ' num2str(times(i)) 's'])
    disp(['  stephen ' num2str(stephentime) 's'])
    disp(['  CVX     ' num2str(t2) 's'])
    disp(['  yCvx =  ' num2str(yCvx')])
    disp(['  y    =  ' num2str(ynic')])
    disp(['  y2   =  ' num2str(ystephen')])
    relErr = norm(ynic - yCvx) / norm(yCvx);
    disp(['  relErr = ' num2str(relErr*100) '%'])
    fprintf('  L(yCvx)  = %.4g\n', L(yCvx))
    fprintf('  L(ystephen)    = %.4g\n', L(ystephen)) 
    fprintf('  L(ynic)    = %.4g\n', L(ynic)) 
end


function x = prox(xbar, d, u)
if isempty(u)
     x = max(xbar, 0);
     return 
end
if all(xbar >= 0)
    x = xbar;
    return 
end
% Checked that this is implemented correctly
b = 1 ./ d;
denom = sqrt(1 + dot(u .* b, u));
v = b .* u / denom;
% Double checked this with the thm
alphas = -xbar ./ (d .* v);
alphas = sort(alphas);
N = length(xbar);
Ls = zeros(N,1);
for n = 1:N
    alpha = alphas(n);
    lambda = xbar + alpha * d .* v;
    mask = lambda > 0;
    Ls(n) = alpha + dot(v, xbar) - dot(v(mask), lambda(mask));
end
[ix, ~] = find(Ls < 0, 1, 'last');
if isempty(ix) % This means that all the constraints are active.
    alphastar = dot(v,xbar);
else
    alpha = alphas(ix);
    lambda = xbar + alpha * d .* v;
    mask = lambda > 0 | (lambda == 0 & d .* v > 0);
    alphastar = (dot(v,xbar) - dot(v(mask), xbar(mask))) / (dot(v(mask), d(mask).*v(mask)) - 1 );
end

% idiot checks
lambdastar = xbar + alphastar * d .* v;
assert(abs(alphastar + dot(v, xbar) - dot(v(mask), lambdastar(mask))) < 1e-12 )
% assert(alphas(ix) < alphastar && (ix == length(d) || alphastar < alphas(ix+1)) ) 

x = max(xbar + alphastar*d.*v, 0);

end % prox
% gcf
% clf
% hold on
% plot(Ns,times)
% xlabel("Ns")
% ylabel("Wall time (s)")