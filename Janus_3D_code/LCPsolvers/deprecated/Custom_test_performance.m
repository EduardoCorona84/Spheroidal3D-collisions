% Copyright, 2011, Kenny Erleben, DIKU.
clear all;
close all;
clc;

max_iter = 10000;
tol_rel  = 1e-5;
tol_abs  = 1e-5;

SizeParameter = 4; % larger than 5 may took ~ 1 hour

%%%%%%%%%%%%%%%%%%%%
%  Fluid Matrix    %
%%%%%%%%%%%%%%%%%%%%

P = []; % Problem size
T = []; % Timings
G = 0;  

for i=1:5
  
  G = G + SizeParameter;
  
  display(['Performance test fluid LCP with G=' num2str(G)])
  
  A      = make_fluid_matrix(G);
  [x, b] = make_lcp(A,0.25);
  x0     = zeros(size(x));
    
  t_start = tic;
  [z, err, iter, flag, conv, m, t0] = multilevel(A, b, x0, 2, max_iter, tol_rel, tol_abs, false);
  t1 = toc(t_start) - t0;
  display(strcat('  Multilevel took ', num2str(t1) ))
  
  t_start = tic;
  fischer_newton(A, b, x0, max_iter, tol_rel, tol_abs, 'random',false);
  t2 = toc(t_start);
  display(strcat('  Fischer took ', num2str(t2) ))
    
  t_start = tic;
  minmap_newton(A, b, x0, max_iter, tol_rel, tol_abs, false);
  t3 = toc(t_start);
  display(strcat('  Min map took ', num2str(t3) ))
  
  t_start = tic;
  pgs(A, b, x0, max_iter, tol_rel, tol_abs, false);
  t4 = toc(t_start);
  display(strcat('  PGS took ', num2str(t4) ))
  
  t_start = tic;
  psor(A, b, x0, 1.4, max_iter, tol_rel, tol_abs, false);
  t5 = toc(t_start);
  display(strcat('  PSOR took ', num2str(t5) ))
  
  t_start = tic;
  BBPGD(A, b, x0, max_iter, tol_rel, tol_abs, false);
  t6 = toc(t_start);
  display(strcat('  BBPGD took ', num2str(t6) ))
  
  t_start = tic;
  APGD(A, b, x0, max_iter, tol_rel, tol_abs, false);
  t7 = toc(t_start);
  display(strcat('  APGD took ', num2str(t7) ))
  
  
  P = [P; length(x)];
  T = [T; [t1, t2, t3, t4, t5, t6, t7] ]; 
end

figure(1)
clf;
set(gca,'FontSize',18);
h1 = semilogy(P, T(:,1), '-','LineWidth',2,'Color',[0.7, 0.1, 0.1]);
grid on;
hold on;
h2 = semilogy(P, T(:,2), ':','LineWidth',2,'Color',[0.1, 0.7, 0.1]);
h3 = semilogy(P, T(:,3), '-.','LineWidth',2,'Color',[0.1, 0.1, 0.7]);
h4 = semilogy(P, T(:,4), '--','LineWidth',2,'Color',[0.7, 0.1, 0.7]);
h5 = semilogy(P, T(:,5), ':','LineWidth',2,'Color',[0.1, 0.7, 0.7]);
h6 = semilogy(P, T(:,6), '-o','LineWidth',2,'Color',[1.0, 0.2, 0.0]);
h7 = semilogy(P, T(:,7), '-o','LineWidth',2,'Color',[0.2, 1.0, 0.2]);
title('Performance Comparison: fluid LCP','FontSize',18);
xlabel('#Variables','FontSize',18);
ylabel('Time [secs]','FontSize',18);
legend([h1, h2, h3, h4, h5, h6, h7], 'Multilevel','Fischer','Min Map', 'PGS', 'PSOR', 'BBPGD', 'APGD');
hold off;
print('-f1','-depsc2','output/performance_fluid');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Random Blocked Matrix    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P = []; % Problem size
T = []; % Timings
G = 0;  

for i=1:5
  
  G = G + SizeParameter;
  
  display(['Performance test random blocked LCP with G=' num2str(G)])
  
  A      = make_blocked_matrix(G,80,0.5);
  [x, b] = make_lcp(A,0.5);
  x0     = zeros(size(x));
    
  t_start = tic;
  [z, err, iter, flag, conv, m, t0] = multilevel(A, b, x0, 2, max_iter, tol_rel, tol_abs, false);
  t1 = toc(t_start) - t0;
  display(strcat('  Multilevel took ', num2str(t1) ))
  
  t_start = tic;
  fischer_newton(A, b, x0, max_iter, tol_rel, tol_abs, 'random',false);
  t2 = toc(t_start);
  display(strcat('  Fischer took ', num2str(t2) ))
    
  t_start = tic;
  minmap_newton(A, b, x0, max_iter, tol_rel, tol_abs, false);
  t3 = toc(t_start);
  display(strcat('  Min map took ', num2str(t3) ))
  
  t_start = tic;
  pgs(A, b, x0, max_iter, tol_rel, tol_abs, false);
  t4 = toc(t_start);
  display(strcat('  PGS took ', num2str(t4) ))
  
  t_start = tic;
  psor(A, b, x0, 1.4, max_iter, tol_rel, tol_abs, false);
  t5 = toc(t_start);
  display(strcat('  PSOR took ', num2str(t5) ))
  
  t_start = tic;
  BBPGD(A, b, x0, max_iter, tol_rel, tol_abs, false);
  t6 = toc(t_start);
  display(strcat('  BBPGD took ', num2str(t6) ))
  
  t_start = tic;
  APGD(A, b, x0, max_iter, tol_rel, tol_abs, false);
  t7 = toc(t_start);
  display(strcat('  APGD took ', num2str(t7) ))
  
  
  P = [P; length(x)];
  T = [T; [t1, t2, t3, t4, t5, t6, t7] ]; 
end

figure(2)
clf;
set(gca,'FontSize',18);
h1 = semilogy(P, T(:,1), '-','LineWidth',2,'Color',[0.7, 0.1, 0.1]);
grid on;
hold on;
h2 = semilogy(P, T(:,2), ':','LineWidth',2,'Color',[0.1, 0.7, 0.1]);
h3 = semilogy(P, T(:,3), '-.','LineWidth',2,'Color',[0.1, 0.1, 0.7]);
h4 = semilogy(P, T(:,4), '--','LineWidth',2,'Color',[0.7, 0.1, 0.7]);
h5 = semilogy(P, T(:,5), ':','LineWidth',2,'Color',[0.1, 0.7, 0.7]);
h6 = semilogy(P, T(:,6), '-o','LineWidth',2,'Color',[1.0, 0.2, 0.0]);
h7 = semilogy(P, T(:,7), '-o','LineWidth',2,'Color',[0.2, 1.0, 0.2]);
title('Performance Comparison: blocked LCP','FontSize',18);
xlabel('#Variables','FontSize',18);
ylabel('Time [secs]','FontSize',18);
legend([h1, h2, h3, h4, h5, h6, h7], 'Multilevel','Fischer','Min Map', 'PGS', 'PSOR', 'BBPGD', 'APGD');
hold off;
print('-f1','-depsc2','output/performance_blocked');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MM-Newton, BBPGD, APGD all fail in this contact matrix unless the matrix
% is small
% MM-Newton stuck in nondescent direction
% BBPGD/APGD stagnation with small step size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P = []; % Problem size
T = []; % Timings
K = 0;  

for i=1:5
  
  K = K + 5;
  display(['Performance test contact LCP with K=' num2str(K)])
  
  A      = make_contact_matrix(K);
  [x, b] = make_lcp(A,0.25);  
  x0     = zeros(size(x));
    
  t_start = tic;
  fischer_newton(A, b, x0, max_iter, tol_rel, tol_abs, 'random',false);
  t1 = toc(t_start);
  display(strcat('  Fischer took ', num2str(t1) ))
  
  t_start = tic;
  minmap_newton(A, b, x0, max_iter, tol_rel, tol_abs, false);
  t2 = toc(t_start);  
  display(strcat('  Min map took ', num2str(t2) ))

  t_start = tic;
  BBPGD(A, b, x0, max_iter, tol_rel, tol_abs, false);
  t3 = toc(t_start);  
  display(strcat('  BBPGD took ', num2str(t3) ))
  
  t_start = tic;
  APGD(A, b, x0, max_iter, tol_rel, tol_abs, false);
  t4 = toc(t_start);  
  display(strcat('  APGD took ', num2str(t4) ))


  
  P = [P; length(x)];
  T = [T; [t1, t2, t3, t4] ]; 
end

figure(3)
clf;
set(gca,'FontSize',18);
h1 = semilogy(P, T(:,1), '-','LineWidth',2,'Color',[0.7, 0.1, 0.1]);
grid on;
hold on;
h2 = semilogy(P, T(:,2), ':','LineWidth',2,'Color',[0.1, 0.7, 0.1]);
h3 = semilogy(P, T(:,3), '-o','LineWidth',2,'Color',[1.0, 0.2, 0.0]);
h4 = semilogy(P, T(:,4), '-o','LineWidth',2,'Color',[0.2, 1.0, 0.2]);
title('Performance Comparison: contact LCP','FontSize',18);
xlabel('#Variables','FontSize',18);
ylabel('Time [secs]','FontSize',18);
legend([h1, h2, h3, h4],'Fischer','Min Map', 'BBPGD', 'APGD');
hold off;
print('-f2','-depsc2','output/performance_contact');