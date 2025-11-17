% In this script we will verify that my implementation of the alternating high and low quasi-Newton (with no inner steps) is the same as Nic's proximal quasi-Newton.

% Set up some sort of LCP
addpath('../../utilities')
addpath('../src')

n = 100;
matrices = construct_test_matrices(n);
A = matrices.linear50.matrix;
global b;
b = randn(n,1);
x0 = zeros(n,1);
clear opts;
opts.A = @(x) A*x;
opts.b = b;
opts.stepSize.init = 'uniform';
opts.stepSize.kappa = 'uniform';
opts.stepSize.eta = 'opt';
opts.storeIts = true;

[x_nic, info_nic] = proxQuasiNewton(@fg, x0, opts);
  
clear opts;
opts.outer = struct();
opts.outer.A = @(x) A*x;
opts.outer.b = b;
opts.outer.stepSize.init = 'uniform';
opts.outer.stepSize.kappa = 'uniform';
opts.outer.stepSize.eta = 'opt';
opts.outer.storeIts = true;
opts.inner.enabled = false; 

[x_mine, info_mine] = alternating_high_and_low(@fg, x0, opts);
disp(info_mine.iterHist - info_nic.iterHist);
disp(norm(x_mine - x_nic));

function [f, grad] = fg(x, Ax)
    global b;
    f = 0.5 * x' * Ax + b' * x;
    grad = Ax + b;
end