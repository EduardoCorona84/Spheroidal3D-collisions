%%
try %#ok<TRYNC>
    rng('default')
end
rng(2);
mfilePath = mfilename('fullpath');
if contains(mfilePath,'LiveEditorEvaluationHelper')
    mfilePath = matlab.desktop.editor.getActiveFilename;
end
[dirname, ~,~] = fileparts(mfilePath);
addpath(genpath(fileparts(dirname)))
fname = 'bimetallic_3x3x3';
load([fname '.mat'], ...
    'A_list', 'b_list');
max_iter = 100;
tol_rel = 1e-6;
tol_abs = 1e-6;
profile = false;
MC = length(A_list);
cvxTime      = zeros(MC,1);
err_bbpgd    = zeros(MC,1);
t_bbpgd      = zeros(MC,1);
iter_bbpgd   = zeros(MC,1);
err_lbfgsb   = zeros(MC,1);
t_lbfgsb     = zeros(MC,1);
iter_lbfgsb  = zeros(MC,1);
err_plbfgs   = zeros(MC,1);
t_plbfgs     = zeros(MC,1);
iter_plbfgs  = zeros(MC,1);
err_zerosr1  = zeros(MC,1);
t_zerosr1    = zeros(MC,1);
iter_zerosr1 = zeros(MC,1);
err_nic      = zeros(MC,1);
t_nic        = zeros(MC,1);
iter_nic     = zeros(MC,1);
err_pq       = zeros(MC,1);
t_pq         = zeros(MC,1);
iter_pq      = zeros(MC,1);
mcGood = [];
for mc = 1:MC
    disp(['mc = ' num2str(mc)])
    Amat = A_list{mc};
    bvec = b_list{mc};
    %% CVX
    N = size(Amat,2);
    x0 = zeros(N,1);
    tic()
    try
        cvx_begin quiet
                variable xRef(N)
                minimize 1/2*dot(xRef, Amat*xRef) + dot(xRef,bvec) 
                subject to 
                0 <= xRef
        cvx_end 
        cvxTime(mc) = toc();
    catch 
        warning('cvx failed probably for negative definiteness of Amat')
        continue
    end
    
    mcGood = [mcGood mc]; %#ok<AGROW>
    nrmXref = norm(xRef);
    errFcn  = @(x) norm( x - xRef )/nrmXref;
    
    %% 'BBPGD'
    tic
    [x_bbpgd, ~ ,iter_bbpgd(mc), ~, ~, ~] = ...
    BBPGD(Amat, bvec, x0, max_iter, tol_rel, tol_abs, profile );    
    t_bbpgd(mc) = toc;
    err_bbpgd(mc) = norm(x_bbpgd - xRef) / norm(xRef);
    %% 'L-BFGS-B'
    tic
    [x_lbfgsb, ~ ,iter_lbfgsb(mc), ~, ~, ~] = ...
    L_BFGS_B(Amat, bvec, x0, max_iter, tol_rel, tol_abs, profile );    
    t_lbfgsb(mc) = toc;
    err_lbfgsb(mc) = norm(x_lbfgsb - xRef) / norm(xRef);
    %% 'P-L-BFGS'
    tic
    [x_plbfgs, ~ ,iter_plbfgs(mc), ~, ~, ~] = ...
    P_L_BFGS(Amat, bvec, x0, max_iter, tol_rel, tol_abs, profile );    
    t_plbfgs(mc) = toc;
    err_plbfgs(mc) = norm(x_plbfgs - xRef) / norm(xRef);
    %% 'zeroSR1'
    tic
    [x_zerosr1, ~, iter_zerosr1(mc), ~, errStruct_stephen, ~] = ...
    ZERO_SR1(Amat, bvec, x0, max_iter, tol_rel, tol_abs, true );   
    t_zerosr1(mc) = toc;
    err_zerosr1(mc) = norm(x_zerosr1 - xRef) / norm(xRef);
    %% 'my zerosr1' 
    fcnGrad = @(x) quadprog(x,Amat,-bvec);
    tic
    opts = struct( ...
        'max_iter',max_iter, ...
        'tol_rel',tol_rel, ...
        'tol_abs',tol_abs, ... 
        'Q', Amat);
    [x_nic, iter_nic(mc), errStruct_nic] = zeroSr1_nic(fcnGrad, x0, opts);
    t_nic(mc) = toc;
    err_nic(mc) = norm(x_nic - xRef) / norm(xRef);
    %% proxQuasiNewton 
    tic
    opts.r = 20;
    [x_pq, iter_pq(mc), errStruct_pq] = proxQuasiNewton(fcnGrad, x0, opts);
    t_pq(mc) = toc;
    err_pq(mc) = norm(x_pq - xRef) / norm(xRef);
    assert(iter_pq(mc) < max_iter);
end
%%
fprintf('Algo     | Rel Err | Time      | Iter\n')
fprintf('CVX      |  N/A    | %.1e s | N/A\n', mean(cvxTime));
fprintf('BB-PGD   |  %.2g  | %.1e s | %.2g\n', ...
    mean(err_bbpgd), mean(t_bbpgd), mean(iter_bbpgd));
fprintf('L-BFGS-B |  %.2g  | %.1e s | %.2g\n', ...
    mean(err_lbfgsb), mean(t_lbfgsb), mean(iter_lbfgsb));
fprintf('P-L-BFGS |  %.2g  | %.1e s | %.2g\n', ...
    mean(err_plbfgs), mean(t_plbfgs), mean(iter_plbfgs));
fprintf('ZEROSR1  |  %.2g  | %.1e s | %.2g\n', ... 
    mean(err_zerosr1), mean(t_zerosr1), mean(iter_zerosr1));
fprintf('nic      |  %.2g  | %.1e s | %.2g\n', ...
    mean(err_nic), mean(t_nic), mean(iter_nic));
fprintf('proxQN   |  %.2g  | %.1e s | %.2g\n', ...
    mean(err_pq), mean(t_pq), mean(iter_pq));
%% 
f = figure();
hold on
edges = [1:1:9 10:5:50, 100];
[N,edges] = histcounts(iter_bbpgd(mcGood),edges);
[N1,edges] = histcounts(iter_lbfgsb(mcGood),edges);
[N2,edges] = histcounts(iter_plbfgs(mcGood),edges);
[N3,edges] = histcounts(iter_zerosr1(mcGood),edges);
[N4,edges] = histcounts(iter_pq(mcGood),edges);
bar([N;N1;N2;N3;N4]');
xlabel('Number of Iterations');
ylabel('Occurance');
tt = join(split(fname, '_'), ' '); 
tt = tt{1};
title(['Comparing LCP solvers for ' tt])
xticks(1:length(edges))
xticklabels([string(edges(1:9)), (string(edges(10:end-1)) + "-" +string(edges(11:end)))])
legend({'BB-PGD', 'L-BFGS-B', 'P-L-BFGS' 'zerosr1', 'proxQuasNewton'})
saveas(f, ['/Users/niru8088/scratch/Spheroidal3D-collisions/docs/' fname '.pdf'])