function xstar = prox(x, opts)
    if isempty(opts.U) && isempty(opts.V)
        % Project with respect to the identity
        xstar = max(0, x);
    elseif size(opts.U,2) + size(opts.V,2) == 1
        % The sign on sigma is counter intuitive, but remember B = B0 + UU' - VV'
        if ~isempty(opts.U)
            sigma = -1;
            w = opts.U; 
        else 
            sigma = 1;
            w = opts.V;
        end
        xstar = prox_rank1(x, opts.h0, w, sigma, opts);
        return  
    else
        xstar = prox_rankr(x, opts.h0, opts.U, opts.V, opts);
    end
end 