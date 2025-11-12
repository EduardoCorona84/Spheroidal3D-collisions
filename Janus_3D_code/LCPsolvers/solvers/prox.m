function xhat = prox(xtilde, opts)
U = opts.prox.U; 
V = opts.prox.V;
if isempty(U) && isempty(V)
    % Project with respect to the identity
    xhat = max(0, xtilde);
elseif size(U,2) + size(V,2) == 1
    % The sign on sigma is counter intuitive, but remember B = B0 + UU' - VV'
    if ~isempty(U)
        sigma = -1;
        w = U; 
    else 
        sigma = 1;
        w = V;
    end
    xhat = prox_rank1(xtilde, opts.prox.H0, w, sigma, opts);
    return  
else
    xhat = prox_rankr(xtilde, opts);
end
end % prox