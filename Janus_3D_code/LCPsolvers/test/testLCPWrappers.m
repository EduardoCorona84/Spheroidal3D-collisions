%%
mfilePath = mfilename('fullpath');
if contains(mfilePath,'LiveEditorEvaluationHelper')
    mfilePath = matlab.desktop.editor.getActiveFilename;
end
[dirname, ~,~] = fileparts(mfilePath);
[dirname, ~] = fileparts(dirname);
addpath(genpath(dirname))
load('LCP_test_case_nic.mat', ...
    'Amat', 'bvec', 'x0', ...
    'max_iter', 'tol_rel', 'tol_abs', ...
    'profile');
%% CVX
N = size(Amat,2);
cvx_begin
        variable xRef(N)
        minimize 1/2*dot(xRef, Amat*xRef) + dot(xRef,bvec) 
        subject to 
        0 <= xRef
cvx_end 
nrmXref = norm(xRef);
errFcn  = @(x) norm( x - xRef )/nrmXref;

%% 'BBPGD'
tic
[x_bbpgd, ~ ,iter_bbpgd, ~, ~, ~] = ...
BBPGD(Amat, bvec, x0, max_iter, tol_rel, tol_abs, profile );    
t_bbpgd = toc;
err_bbpgd = norm(x_bbpgd - xRef) / norm(xRef);
%% 'L-BFGS-B'
tic
[x_lbfgsb, ~ ,iter_lbfgsb, ~, ~, ~] = ...
L_BFGS_B(Amat, bvec, x0, max_iter, tol_rel, tol_abs, profile );    
t_lbfgsb = toc;
err_lbfgsb = norm(x_lbfgsb - xRef) / norm(xRef);
%% 'P-L-BFGS'
tic
[x_plbfgs, ~ ,iter_plbfgs, ~, ~, ~] = ...
P_L_BFGS(Amat, bvec, x0, max_iter, tol_rel, tol_abs, profile );    
t_plbfgs = toc;
err_plbfgs = norm(x_plbfgs - xRef) / norm(xRef);
%% 'zeroSR1'
tic
[x_zerosr1, ~, iter_zerosr1, ~, ~, ~] = ...
ZERO_SR1(Amat, bvec, x0, max_iter, tol_rel, tol_abs, profile );   
t_zerosr1 = toc;
err_zerosr1 = norm(x_zerosr1 - xRef) / norm(xRef);
%%
fprintf('Algo     | Rel Err | Time | Iter\n')
fprintf('BBGPD    | %.2g | %.1g s | %d\n', err_bbpgd, t_bbpgd, iter_bbpgd);
fprintf('L-BFGS-B | %.2g | %.1g s | %d\n', err_lbfgsb, t_lbfgsb, iter_lbfgsb);
fprintf('P-L-BFGS | %.2g | %.1g s | %d\n', err_plbfgs, t_plbfgs, iter_plbfgs);
fprintf('ZEROSR1  | %.2g | %.1g s | %d\n', err_zerosr1, t_zerosr1, iter_zerosr1);

%%
Amatvec = @(x) Amat*x; 
%% 'BBPGD'
tic
[x_bbpgd, ~ ,iter_bbpgd, ~, ~, ~] = ...
BBPGD_matfree(Amatvec, bvec, x0, max_iter, tol_rel, tol_abs, profile );    
t_bbpgd = toc;
err_bbpgd = norm(x_bbpgd - xRef) / norm(xRef);
%% 'L-BFGS-B'
tic
[x_lbfgsb, ~ ,iter_lbfgsb, ~, ~, ~] = ...
L_BFGS_B_matfree(Amatvec, bvec, x0, max_iter, tol_rel, tol_abs, profile );    
t_lbfgsb = toc;
err_lbfgsb = norm(x_lbfgsb - xRef) / norm(xRef);
%% 'P-L-BFGS'
tic
[x_plbfgs, ~ ,iter_plbfgs, ~, ~, ~] = ...
P_L_BFGS_matfree(Amatvec, bvec, x0, max_iter, tol_rel, tol_abs, profile );    
t_plbfgs = toc;
err_plbfgs = norm(x_plbfgs - xRef) / norm(xRef);
%% 'zeroSR1'
tic
[x_zerosr1, ~, iter_zerosr1, ~, ~, ~] = ...
ZERO_SR1_matfree(Amatvec, bvec, x0, max_iter, tol_rel, tol_abs, profile );   
t_zerosr1 = toc;
err_zerosr1 = norm(x_zerosr1 - xRef) / norm(xRef);
%%
disp('============ MAT FREE ===============')
fprintf('Algo     | Rel Err | Time | Iter\n')
fprintf('BBGPD    | %.2g | %.1g s | %d\n', err_bbpgd, t_bbpgd, iter_bbpgd);
fprintf('L-BFGS-B | %.2g | %.1g s | %d\n', err_lbfgsb, t_lbfgsb, iter_lbfgsb);
fprintf('P-L-BFGS | %.2g | %.1g s | %d\n', err_plbfgs, t_plbfgs, iter_plbfgs);
fprintf('ZEROSR1  | %.2g | %.1g s | %d\n', err_zerosr1, t_zerosr1, iter_zerosr1);