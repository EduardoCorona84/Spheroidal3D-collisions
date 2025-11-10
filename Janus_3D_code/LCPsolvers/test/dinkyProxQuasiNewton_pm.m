%% Set up problem 
rng(1);
n = 100; % size of the original problem 
r1 = 3; % rank of the positive update
r2 = 0; % rank of the neg update
r = r1 + r2;
d = 100;
MC = 100;
runTime = zeros(MC,1);
cvxTime = zeros(MC,1);
relErr = zeros(MC,1);
for mc = 1:MC
    disp(['mc ' num2str(mc)])
    U = rand(n,r1);
    V = .5*rand(n,r2);
    B = eye(n)*d + U*U' - V*V';
    y = randn(n,1);
    % custom solver
    tic 
    xstar = prox(y, d, U, V);
    runTime(mc) = toc;
    % cvx
    tic
    f = @(x)1/2*dot(x-y, B*(x-y)) ;
    cvx_begin quiet
            variable xRef(n)
            minimize f(xRef)
            subject to 
            0 <= xRef
    cvx_end 
    cvxTime(mc) = toc;
    relErr(mc) = norm(xRef - xstar) / norm(xRef);
end
%% compare
disp(['Problem Size: n = ' num2str(n) ', r = ' num2str(r)])
disp('Custom Solver')
disp(['  mean run time ' num2str(mean(runTime)) 's, std ' num2str(std(runTime)) 's'])
disp('  xstar')
disp(['    ' num2str(xstar(xstar>0)')])
disp('CVX')
disp(['  mean run time ' num2str(mean(cvxTime)) 's, std ', num2str(std(cvxTime)) 's'])
disp('  xRef')
disp(['    ' num2str(xRef(xstar>0)')])
disp('Relative Error')
disp([ '  mean ' num2str(mean(relErr)) ' std ' num2str(std(relErr))])