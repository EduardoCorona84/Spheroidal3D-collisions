function [x_k, info] = back_track(fg, x_k, x_km1, kappa, p, opts, info)
    %Perform back tracking line search.
    %For now we only will do a non monotone line search.
    %basically ripped from FASTA paper
    if opts.acceleration.backtrack
        f_max = max(opts.acceleration.f_vals);
        [f_k, ~] = fg(x_k);
        while f_k - 1e-12 > f_max + -(x_k - x_km1)'*p + (1/(2*kappa))*norm(x_k - x_km1)^2
            kappa = kappa / 2;
            x_k = x_km1 + kappa * p;
            x_k = max(0, x_k);
            [f_k, ~] = fg(x_k);
        end
        opts.acceleration.f_vals = [opts.acceleration.f_vals(2:end); f_k];
        return
    
    else
        return
    end