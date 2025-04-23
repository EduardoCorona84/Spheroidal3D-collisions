function [x, flag] = line_search(out, fgFcn, options)
% Stephen adding fx_Ax and A inputs, and Ax output
%%
%%  Armijo along projection arc.
%%

step = 1;
x = out.x + step * out.srch;
flag = -1;

for i = 1 : options.max_func_evals
    x(x < 0) = 0;
    delx = out.x - x;
    fc = fgFcn(x);

    if out.obj - fc >= options.sigma * out.grad' * delx
        flag = 1;
        return;
    end

    step = step * options.beta;
    x = out.x + step * out.srch;
end

x = out.x;  % we should never get to this line
