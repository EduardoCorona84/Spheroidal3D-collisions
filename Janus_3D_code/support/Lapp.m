
function y = Lapp(A,x, verbose)
if ~exist('verbose','var') || isempty(verbose)
    verbose = false;
end
if verbose
    tic
end
if isnumeric(A)
    y=A*x; 
else
    y=real(A(x)); 
end
if verbose
    dt = toc; 
    fprintf('\n -- mvp time %.4g sec', dt)
end
end
