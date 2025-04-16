plotFlag = true;
Ns = arrayfun(@(e) round(10^e), 1:0.5:4);
times = zeros(length(Ns));
for i = 1:1%length(Ns)
    N = Ns(i);
    x = 100*randn(N,1);
    d = abs(100*randn(N,1));
    u = 100*randn(N,1);
    
    tic
    y = prox_rank1_nn(x, d, u, plotFlag);
    times(i) = toc();

    res = x-y;
    L = @(z) 1/2*sum((x-z).*d.*(x-z)) ...
        + dot(u,(x-z))^2;
    % tic
    % cvx_begin quiet
    %     variable yCvx(N)
    %     minimize L(yCvx)
    %     subject to 
    %         0 <= yCvx
    % 
    % cvx_end
    % t2 = toc;
    disp(['N = ' num2str(N)])
    disp(['  My Algo ' num2str(times(i)) 's'])
    % disp(['  CVX ' num2str(t2) 's'])
    % disp(['  yCvx = ' num2str(yCvx')])
    % disp(['  y = ' num2str(y')])
    % relErr = norm(y - yCvx) / norm(yCvx);
    % disp(['  relErr = ' num2str(relErr*100) '%'])
    % fprintf('  L(yCvx) = %.4g\n', L(yCvx))
    fprintf('  L(ystar) = %.4g\n', L(y)) 
end

% gcf
% clf
% hold on
% plot(Ns,times)
% xlabel("Ns")
% ylabel("Wall time (s)")