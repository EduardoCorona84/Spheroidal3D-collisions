function plot_singular_values(mtx, tol)
    [~, S, ~] = svd(mtx);
    s_vals = diag(S);

    rank_no_tol = rank(mtx);
    rank_tol = rank(mtx, tol);
    fprintf("Rank (no given tolerance): %d -- nullspace dimension is %d\n", rank(mtx), size(mtx,1) - rank_no_tol);
    fprintf("Rank (tolerance of %f): %d -- nullspace dimension is %d\n", tol, rank_tol, size(mtx,1) - rank_tol);
    
    figure;
    semilogy(s_vals, 'o-');
    title('Singular Values of Matrix');
    xlabel('Index');
    ylabel('Singular Value');
    grid on;
end