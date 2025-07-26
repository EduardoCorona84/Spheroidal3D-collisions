function out = pqn_general(fgFcn, x0, opt)
%
% function out = pqn_general(fgFcn, x0, opt)
% Solve a convex problem on the nonnegative orthant,
% min f, s.t. x \ge 0
%
% fgFcn --
% x0    -- Starting vector (useful for warm-starts).
% OPT   -- This structure contains important opt that control how the
%          optimization procedure runs. To obtain a default structure 
%          the user can use 'opt = solopt'. Use 'help solopt' to get a 
%          description of what the individual opt mean.
%
% OUT   -- contains the solution and other information about the 
%          optimization or info indicating whether the method succeeded or 
%          failed.
%
%
% Version 1.2 (c) 2009  Dongmin Kim  and Suvrit Sra
% Modified April 2025 by Nic Rummel 

out = [];
if strcmpi(opt.algo, 'PLB') 
    out = plb_general(fgFcn, x0, opt);
    return 
end
    
error([opt.algo ' is not extended to the quadratic case'])

    
