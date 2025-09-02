function [converged, opts, info] = check_convergence(k, f_k, f_km1, x_k, x_km1, info, opts)
    %in this function, we have to do 3 things
    %1. Compare the recent high and low fidelity iterates and make a determination regarding if we have sufficient decrease (under some criterion). If we don't have sufficient decrease we want to adapt the high/low fidelity split in some way
    %this is under opts.outer.adaptive
    %2. If we are in the high fidelity only regime (or some to be developed equivalent for alternating high low), we want to check if we have converged using various criterion (ideas rel, abs, kkt, with rel and abs (still want to characterize this))
    %3.Update the info struct with the information used to make these decisions

    if k == 0
        %for now we will do nothing. If KKT/merit function is a good measure, we could use this here.
    else
        if k >= opts.outer.max_iter
            converged = true;
        else
            if opts.inner.enabled == false
            %we are just checking convergence and updating info
                abs_error = norm(x_k - x_km1);
                if abs_error < opts.outer.solver_opts.tol_abs
                    converged = true;

                elseif abs_error/norm(x_k) < opts.outer.solver_opts.tol_rel
                    converged = true;
                end

            elseif opts.inner.enabled == true
                if f_k < f_km1
                    %do nothing
                else
                    switch(lower(opts.outer.adaptive))
                        case 'high'
                            opts.inner.enabled = false;

                        case 'halving'
                            if opts.inner.solver_opts.max_iter == 1
                                opts.inner.enabled = false;
                            else
                                opts.inner.solver_opts.max_iter = floor(opts.inner.solver_opts.max_iter / 2);
                            end

                    end
                end
            end
        end
        
        if opts.outer.storeIts == true
            info.outer.iterHist(k+1, :) = x_k;
        end

    end

    if converged == true
        info.outer.iter = k;
    end

end







    %this is nic's check converence. I have some questions on this, so just gonna make a barebones version for now
    %{
    persistent old_kkt
    if k == 0
        old_kkt = [];
    end
    % Check for convergence
    converged = false;
    phi = min(grad_k,x_k);
    kkt = 0.5*dot(phi, phi);
    if k >= opts.max_iter
        info.flag =  8;
        converged = true;
    % elseif eta < 10*eps
    %     % if step direction gets too small give up
    %     info.flag =  5;
    %     converged = true;
    else
        % kkt conditions / LCP being satisified is equivalent to
        % the grad_k(i) = 0 \perp x_k(i) = 0
        if ~isempty(old_kkt) && (abs(kkt - old_kkt) / abs(kkt)) < opts.tol_rel
            % Relative stopping criteria
            info.flag = 3;
            converged = true;
        elseif kkt < opts.tol_abs
            % Absolute stopping criteria
            info.flag = 4;
            converged = true;
        end
        old_kkt = kkt;
    end

    % Keep track of error history
    if ~isempty(opts.errFcn)
        if isa(opts.errFcn,'function_handle')
            info.errHist(k+1) = opts.errFcn(x_k);
        elseif iscell(opts.errFcn)
            for i = 1:numel(opts.errFcn)
                fcn = opts.errFcn{i};
                info.errHist(k+1,i) = fcn(x_k);
            end
        end
    end

    if opts.storeIts 
        info.iterHist(k+1,:) = x_k;
    end
    % Fill up info if converged
    if converged 
        if k == 0 
            info.flag = 9;
        end
        info.iter = k;
        info.f = f_k;
        info.kkt = kkt;
        info.msg = flag2msg(info.flag);
        if ~isempty(opts.errFcn)
            info.errHist = info.errHist(1:k+1,:);  
        end
        if opts.storeIts
            info.iterHist = info.iterHist(1:k+1,:);
        end
    end
end % checkConvergence

function msg = flag2msg(flag)
    assert(1 <= flag && flag <= 9)
    % Just a list of human readable text strings to convert the flag return
    % code into something readable by writing msg(flag) onto the screen.
    msgs = {...
        'preprocessing';  % info.flag =  1
        'iterating';      % info.flag =  2
        'relative';       % info.flag =  3
        'absolute';       % info.flag =  4
        'stagnation';     % info.flag =  5
        'local minima';   % info.flag =  6
        'nondescent';     % info.flag =  7
        'maxlimit';       % info.flag =  8
        'x0 is sufficient' % info.flag =  9
    };
    msg = msgs{flag};
    
end % flag2msg
%}