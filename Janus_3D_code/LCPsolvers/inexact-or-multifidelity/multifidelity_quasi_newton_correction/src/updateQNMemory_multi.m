function [r, h0, rho, S, Y, opts] = updateQNMemory_multi(k, s, y, opts)
    %this had used persistent variables, but I personally prefer storing these things in opts and then updating
    %persistent rho_ S_ Y_ 

    n = numel(s);
    r = opts.r;
    if k==0 || any([~isfield(opts, 'rho_'), ~isfield(opts, 'S_'), ~isfield(opts, 'Y_')])
        opts.rho_ = zeros(r, 1);
        opts.S_ = zeros(n, r);
        opts.Y_ = zeros(n, r);
    end
    r = min(r, k);
    if r <= 0
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
    % Loose memory if the iterate is larger than memory
    if k > r
        opts.rho_(1:r-1) = opts.rho_(2:r);
        opts.S_(:,1:r-1) = opts.S_(:,2:r);
        opts.Y_(:,1:r-1) = opts.Y_(:,2:r);
        opts.rho_(r) = 0;
        opts.S_(:,r) = 0;
        opts.Y_(:,r) = 0;
    end
    % Set r to the current amount of memory, handeling the case previous
    % hessian updates have been skipped due to curvature checks failing
    r_ = find(arrayfun(@(i) all(opts.S_(:,i) == 0), 1:r),1,'first');
    if ~isempty(r_)
        r = min(r, r_);
    end
    % Curvature check to see if we should use this update to for hessian, %
    %TO DO: Either Here or get_correction_SR1, we need to add additional logic for skipping sr1 updates using condition 6.26 in nocedal and wright.
    rho_candidate = 1/ dot(s,y);
    if 1/rho_candidate >= 1e-8
        opts.S_(:, r) = s;
        opts.Y_(:, r) = y;
        opts.rho_(r) = rho_candidate;
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
    rho = opts.rho_(1:r);
    S = opts.S_(:, 1:r);
    Y = opts.Y_(:, 1:r);
    % Set h0
    tau_bb2 = 1 / rho(r) / dot(Y(:, r),Y(:, r)); % dot(s_k,y_k) / norm(y_k,2)^2
    tau_bb2 = min(max(tau_bb2, opts.tau_min), opts.tau_max);
    if tau_bb2 == opts.tau_min
        warning('Convexity of cost function is stagnating'); 
    end
    h0 = opts.gamma * tau_bb2;
end % updateQNMemory


