function solvePreconditionedProblem(A, b, M, opts)
%% We are given a proconditioner so that M A M \approx I
% We perform a change of variables y = M^{-1} x
B = @(x) M * A(M*x);
c = M \ b; 
% Solve \min_{-M y \leq 0} 1/2 y^T B y + y^Tc
% Unfortunate this no longer has a separable prox for the inequality 
% So we solve the dual problem 



end