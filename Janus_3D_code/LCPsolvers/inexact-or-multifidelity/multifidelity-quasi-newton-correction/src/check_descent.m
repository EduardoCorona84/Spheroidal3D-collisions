function [error, descent] = check_descent(p, opts)
    %This function computes an upper bound on the error in the low fidelity evaluation and checks the descent condition
    

    %We want to compute a bound on the error of (A - A_hat)*p
    %Recall that we have that A_hat p is exact if p is in the span of the columns of S
    %Thus we we want to find the norm of the projection of p onto the orthogonal complement of the span of S
    %Suppose we have some estimate of the norm of A - A_hat, call abs_err
    %Thus the error would be bounded by abs_err * || P_{S_perp} p ||_2

    
