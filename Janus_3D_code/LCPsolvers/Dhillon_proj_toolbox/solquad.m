function out = solnls(A, b, x0, opt)
%
% function out = solnls(A, b, x0, opt)
% Solve a nonnegative least squares problem,
% min 1/2x'Ax - x'b, s.t. x \ge 0
%
% x0 -- Starting vector (useful for warm-starts).
%
% OPT -- This structure contains important opt that control how the
% optimization procedure runs. To obtain a default structure the user can
% use 'opt = solopt'. Use 'help solopt' to get a description of
% what the individual opt mean.
%
% OUT contains the solution and other information about the optimization or
% info indicating whether the method succeeded or failed.
%
% See also: solopt
%
% Version 1.2 (c) 2009  Dongmin Kim  and Suvrit Sra
% 
out = [];
if strcmpi(opt.algo, 'PLB') 
    out = plb_quad(A, b, x0, opt);
    return 
end
    
error([opt.algo ' is not extended to the quadratic case'])

    
