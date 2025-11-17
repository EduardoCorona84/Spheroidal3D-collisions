function [r, h0, rho, S, Y] = updateQNMemory(k, s, y, opts)
    persistent rho_ S_ Y_ 

    n = opts.n;
    r = opts.r;

    %There is an absolute ton of duplicated logic here, but I am just trying to get this to work for now and then I can refactor later.
    %Should add this to check opts, we will only use memory assuming the same r for all problems

    %flag to return the memory without updating it
    if k == -2
        r_ = find(arrayfun(@(i) all(S_(:,i) == 0), 1:r),1,'first');
        if ~isempty(r_)
            r = min(r, r_);
        else
            %do nothing
        end
        r = r - 1;
        h0 = [];
        rho = rho_;
        S = S_;
        Y = Y_;
        % Set h0
        return
    elseif k == 1
        % Initialize memory
        if opts.qn_memory.warm == true
            rho_ = opts.qn_memory.rho;
            S_ = opts.qn_memory.S;
            Y_ = opts.qn_memory.Y;
        else
            rho_ = zeros(r, 1);
            S_ = zeros(n, r);
            Y_ = zeros(n, r);
        end
        %find the first unused slot
    else
        %do nothing
    end
    %find the current memory
    r_ = find(arrayfun(@(i) all(S_(:,i) == 0), 1:r),1,'first');
    %if its empty, we have full memory, find the current amount of memory
    if ~isempty(r_)
        r = r_; %slot we will write to
    else
        r = r; %slot we will write to
    end


    if r <= 1 && k == 1
        % Initial step we do not compute any curvature information because only
        % information about x0 is known, thus no secant equaion could be
        % satisfied
        r = 0;
        h0 = 1; 
        rho = [];
        S = [];
        Y = [];
        return 
    end

    %now we update the memory
    if k == 1
        %do nothing, we don't have a secant equation for this problem yet. 
        %I suppose we could try this by passing the last evaluation of the previous problem?
    else
        %if memory is full, we loose the oldest
        if isempty(r_)
            rho_(1:r-1) = rho_(2:r);
            S_(:,1:r-1) = S_(:,2:r);
            Y_(:,1:r-1) = Y_(:,2:r);
            rho_(r) = 0;
            S_(:,r) = 0;
            Y_(:,r) = 0;
        end
        %now we check if we can add the new info
        rho_candidate = 1/ dot(s,y);
        if 1/rho_candidate >= 1e-8
            S_(:, r) = s;
            Y_(:, r) = y;
            rho_(r) = rho_candidate;
        elseif r == 1
            % No memory is available and current udpate is bad so do the gradient
            % descent step
            r = 0;
            h0 = 1; 
            rho = [];
            S = [];
            Y = [];
            return 
        else 
            % if the curvature check fails use the previous information
            r = r - 1;
        end

        
    end

    rho = rho_(1:r);
    S = S_(:, 1:r);
    Y = Y_(:, 1:r);
    % Set h0
    tau_bb2 = 1 / rho(r) / dot(Y(:, r),Y(:, r)); % dot(s_k,y_k) / norm(y_k,2)^2
    tau_bb2 = min(max(tau_bb2, opts.tau_min), opts.tau_max);
    if tau_bb2 == opts.tau_min
        warning('Convexity of cost function is stagnating'); 
    end
    h0 = opts.gamma * tau_bb2;

end


