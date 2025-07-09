function out = plb_general(fgFcn, x0, options)
% P_LBFGS -- This function solves the following optimization problem
%  min (1/2)x^TAX + x^Tb, subject to x \ge 0
% 
% based on projected L-BFGS method in 'Tackling Box-Constrained
% Optimization via a New Projected Quasi-Newton Approach' by Dongmin et al
% 
% Usage:
%       out = plb_quad(A, b, x0, options)
% IN
%  A  -- matrix or mat-vec operator in the objective fun
%  b  -- vector in the objective fun
%  x0 -- Starting vector (useful for warm-starts) -- *can* be zero.
%  options -- options for optimization (see solopt)
% OUT 
%  out -- contains the solution and other information about the optimization 
%         or info indicating whether the method succeeded or failed.
%
% Version 1.2 (c) 2009  Dongmin Kim  and Suvrit Sra
% Extended from code from Version 1.2 (c) 2009  Dongmin Kim  and Suvrit Sra
% Extended by Stephen Becker for comparison in ~2014
% Extend by April 2024 Nic Rummel for quadratic program instead of nnlsq
%% ------------------------------------------------------
%  INITIALIZATION
%  ------------------------------------------------------
out = getout();
out.time = 0;
out.algo = 'PLB_quad';
out.start_time = clock;
out.status = 'Failure';
h = 0;

%% -------------------------------------------------------
%  Prime the pump
%  -------------------------------------------------------

if (options.verbose)
    if (options.asgui) h=waitbar(0, 'Running PLB...'); 
    else fprintf('PLB:\n'); 
    end
end

%% ------------------------------------------------------
%  OTHER INITIALIZATION
%  ------------------------------------------------------
rho = ones(options.maxmem, 1);
alp = rho;
delx = zeros(length(x0), options.maxmem);
delg = delx;
last = 1;

out.x = x0;
out.oldx = x0;
[out.obj,out.grad]  = fgFcn(x0);
out.oldobj = out.obj;
out.oldgrad = out.grad;
out.srch = -out.grad;

[out.x, ~] = line_search(out, fgFcn, options);
[out.obj, out.grad]  = fgFcn(out.x);
out.iter = 1;
%
if isfield(options, 'errFcn')
    if ishandle(options.errFcn)
        out.errHist = zeros(options.maxit+1,1);
        out.errHist(1) = options.errFcn(x0);
        out.errHist(2) = options.errFcn(x0);
    elseif iscell(options.errFcn)
        out.errHist = zeros(options.maxit+1,numel(options.errFcn));
        for i = 1:numel(options.errFcn)
            fcn = options.errFcn{i};
            out.errHist(1,i) = fcn(x0);
            out.errHist(2,i) = fcn(out.x);
        end
    else 
        assert(false, 'errFcn must be empty, a function handle or a cell array full of function handles');
    end
end
%% -----------------------------------------------------
%  The main iterative loop
%  -----------------------------------------------------
while true
    out.iter = out.iter + 1;
    show_status(out.iter, options, h);

    gp = find(out.x == 0 & out.grad > 0);

    % Compute L-BFGS.
    delx(:, last) = out.x - out.oldx;
    delg(:, last) = out.grad - out.oldgrad;
    delx(gp, :) = 0;
    delg(gp, :) = 0;

    out.oldx = out.x;
    out.oldgrad = out.grad;
    out.oldobj = out.obj;

    out.grad(gp) = 0;
    out.srch = out.grad;
    rho(last) = 1 / (delx(:, last)' * delg(:, last));
    pt = last;

    for i = 1 : min(out.iter, options.maxmem)
        alp(pt) = rho(pt) * delx(:, pt)' * out.srch;
        out.srch = out.srch - alp(pt) * delg(:, pt);
        pt = options.maxmem - mod(-pt + 1, options.maxmem);
    end

    out.srch = 1 / rho(last) / (delg(:, last)' * delg(:, last)) * out.srch;
    
    for i = 1 : min(out.iter, options.maxmem)
        pt = mod(pt, options.maxmem) + 1;
        b = rho(pt) * delg(:, pt)' * out.srch;
        out.srch = out.srch + (alp(pt) - b) * delx(:, pt);
    end

    last = mod(last, options.maxmem) + 1;

    out.srch = -out.srch;
    out.srch(gp) = 0;
    % SRB's linesearch version can make use of the mat vec
    [out.x ,~] = line_search(out, fgFcn, options);
    [out.obj, out.grad] = fgFcn(out.x);
    % Nic adding:
    if isfield(options, 'errFcn')
        if ishandle(options.errFcn)
            out.errHist(out.iter+1) = options.errFcn(out.x);
        elseif iscell(options.errFcn)
            for i = 1:numel(options.errFcn)
                fcn = options.errFcn{i};
                out.errHist(out.iter+1,i) = fcn(out.x);
            end
        end
    end
    % termination
    if flag < 0
        term_reason = 6;
    else
        term_reason = check_termination(options, out);
    end
    if (term_reason > 0)
        break;
    end
%%

end % of while

%% ------------------------------------------------
%  Final statistics and wrap up
%  ------------------------------------------------
out.time = etime(clock, out.start_time);
out.status = 'Success';
out.term_reason = set_term_reason(term_reason);

if (options.verbose)
    if (options.asgui) delete(h); else fprintf('Done\n'); end
end

if isfield(options, 'errFcn')
    out.errHist= out.errHist(1:out.iter+1,:);
end


