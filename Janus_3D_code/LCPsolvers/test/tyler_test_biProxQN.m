%In this script we are going to test the alternating high and low fidelity quasi-Newton method on some data.

%In this script, we will plot the multifidelity quasi-Newton correction for different low fidelities for a single example problem.

addpath('../solvers/');
addpath('../../goodData/');

% Load Data
file_name = 'amphi.special.n_3.p_8.cDist_2.3.allMats.mat';
load(file_name);

%find reference high-fidelity solution
high_fidelity_mat = A{32, length(ps), 1};
b_vec = b{32};

% compute condition number (use condest(high_fidelity_mat) for large/sparse)
cond_high = cond(high_fidelity_mat);

clear opts;

error_function = @(x) dot(min(x, high_fidelity_mat*x + b_vec), min(x, high_fidelity_mat*x + b_vec));
opts.solver = 'proxQuasiNewton';
opts.storeIts = true; 
opts.errFcn = error_function;
opts.A = @(x) high_fidelity_mat*x;
opts.b = b_vec;

fg_high = @(x, Ax) quadraticLoss(x, opts.A, b_vec, Ax);
[x_high, info_high] = proxQuasiNewton(fg_high, zeros(size(b_vec)), opts);



%now find the low fidelities
errors_fixed = cell(1, length(ps));
for p = 1:length(ps)
    low_fidelity_mat = A{32, p, 4};
    clear opts;
    opts.solver = 'bifi';
    opts = defaultLCPOpts(opts, b_vec);
    opts.errFcn = error_function;
    opts.A = @(x) high_fidelity_mat*x;
    opts.b = b_vec;
    opts.low.A = @(x) low_fidelity_mat*x; 
    opts.low.b = b_vec;
    global fixBifi 
    fixBifi = true;
    [x_low, info_low] = bifidelityProxQuasiNewton(fg_high, zeros(size(b_vec)), opts);
    errors_fixed{p} = info_low.errHist;
end

%create figures for different p's. 

%now find the low fidelities
errors = cell(1, length(ps));
for p = 1:length(ps)
    low_fidelity_mat = A{32, p, 4};
    clear opts;
    opts.solver = 'bifi';
    opts = defaultLCPOpts(opts, b_vec);
    opts.errFcn = error_function;
    opts.A = @(x) high_fidelity_mat*x;
    opts.b = b_vec;
    opts.low.A = @(x) low_fidelity_mat*x; 
    opts.low.b = b_vec;
    global fixBifi 
    fixBifi = false;
    [x_low, info_low] = bifidelityProxQuasiNewton(fg_high, zeros(size(b_vec)), opts);
    errors{p} = info_low.errHist;
end

%create figures for different p's. 

cond_high = cond(high_fidelity_mat);
figure;
semilogy(info_high.errHist, '-o', 'DisplayName', sprintf('High Only p=%d, tol=%.1e', ps(length(ps)), tols(4)));
hold on;
for p = 1:length(ps)
    semilogy(errors_fixed{p}, '--o', 'DisplayName', sprintf('Multi Fixed p=%d, tol=%.1e', ps(p), tols(4)));
    semilogy(errors{p}, '--o', 'DisplayName', sprintf('Multi p=%d, tol=%.1e', ps(p), tols(4)));
end
xlabel('Outer Iteration');
ylabel('KKT Error');
title(sprintf('Convergence of Multi-Fidelity proximal Quasi-Newton (cond(A)=%.2e)', cond_high));
legend('show', 'Location', 'best');
saveas(gcf, 'bifidelityProxQN_varyingLowFidelity.svg');
hold off;
