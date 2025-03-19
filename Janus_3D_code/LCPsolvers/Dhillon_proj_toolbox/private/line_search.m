function [x, flag, step, Ax] = line_search(out, fx, gfx, options, fx_Ax, A)
% Stephen adding fx_Ax and A inputs, and Ax output
%%
%%  Armijo along projection arc.
%%

if nargin < 6 || isempty( fx_Ax ) || isempty(A)
    DO_AX = false;
    Ax    = [];
else
    DO_AX = true;
    Ax     = [];
end

    step = 1;
    x = out.x + step * out.srch;
    flag = -1;
    
    for i = 1 : options.max_func_evals
        x(x < 0) = 0;
        delx = out.x - x;
        if DO_AX
            Ax = A*x;
            fc = fx_Ax(Ax);
        else
            fc = fx(x);
        end
        
        if out.obj - fc >= options.sigma * out.grad' * delx
            flag = 1;
            return;
        end
        
        step = step * options.beta;
        x = out.x + step * out.srch;
    end

    x = out.x;  % we should never get to this line
