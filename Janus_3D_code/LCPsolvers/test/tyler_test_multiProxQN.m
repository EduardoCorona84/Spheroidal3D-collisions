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

clear opts;
opts.solver = 'proxquasinewton';
opts.storeIts = true;

error_function = @(x) (1/2)*dot(min(x, high_fidelity_mat*x + b_vec), min(x, high_fidelity_mat*x + b_vec));
opts.errFcn = error_function;
opts.A = @(x) high_fidelity_mat*x;
opts.b = b_vec;

fg_high = @(x, Ax) quadraticLoss(x, opts.A, b_vec, Ax);
[x_high, info_high] = proxQuasiNewton(fg_high, zeros(size(b_vec)), opts);



%now find the low fidelities
errors = cell(length(tols), length(ps));
for row = 1:length(tols)
    for col = 1:length(ps)
        low_fidelity_mat = A{32, col, row};
        clear opts;
        opts.high.solver = 'proxquasinewton';
        opts.high.errFcn = error_function;
        opts.high.storeIts = true;
        opts.high.A = @(x) high_fidelity_mat*x;
        opts.high.b = b_vec;
        opts.high.max_iter = 8;
        opts.sub.solver = 'proxquasinewton';  
        opts.low.A = @(x) low_fidelity_mat*x;
        opts.low.b = b_vec;
        opts.low.max_iter = 4;
        [x_low, info_low] = multifidelityProxQuasiNewton(fg_high, zeros(size(b_vec)), opts);

        errors{row, col} = info_low.errHist;

    end

end

%create figures for different p's. 
for col = 1:length(ps)
    figure;
    semilogy(info_high.errHist, '-o', 'DisplayName', 'High');
    hold on;
    for row = 1:length(tols)
        semilogy(errors{row, col}, '-o', 'DisplayName', sprintf('p=%d, tol=%.1e', ps(col), tols(row)));
    end
    xlabel('Outer Iteration');
    ylabel('KKT Error');
    title(sprintf('No Correction Convergence for p=%d', ps(col)));
    legend('show', 'Location', 'best');
    saveas(gcf, sprintf('no_correction_convergence_p_%d.fig', ps(col)));
    hold off;
end