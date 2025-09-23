function opts = update_correction_memory(k, s_k, y_k, A_low_s_k, opts)
    %this function updates the memory storage for the limited memory/compact representations

    %if we are going to store the updates, we only need to store the S_k and (Y_k - B_0*S_k) matrices. In this case, we let Y = y_k - B_0*s_k

    %we also need to consider if we are storing a correction along b for the entire optimization problem
    if opts.curr_mem == 0
        %first update, just add the first column
        opts.S(:, 1) = s_k;
        %now we check which update we are using, and will store the appropriate vectors in the appropriate matrix
        switch lower(opts.update)
            case 'sr1'
                opts.utility_matrix(:, 1) = y_k - A_low_s_k; 
            case {'bfgs', 'dfp'}
                opts.utility_matrix(:, 1) = A_low_s_k;
                opts.Y(:, 1) = y_k;
        end
        opts.curr_mem = 1;
    elseif opts.curr_mem < opts.r
        %not full memory yet, just add the next column
        opts.S(:, opts.curr_mem + 1) = s_k;
        switch lower(opts.update)
            case 'sr1'
                opts.utility_matrix(:, opts.curr_mem + 1) = y_k - A_low_s_k; 
            case {'bfgs', 'dfp'}
                opts.utility_matrix(:, opts.curr_mem + 1) = A_low_s_k;
                opts.Y(:, opts.curr_mem + 1) = y_k;
                
        end
        opts.curr_mem = opts.curr_mem + 1;
    else
        %full memory, need to shift everything over and add the new column at the end
        if opts.persistent_first_update == true 
            initial_index = 2;
        else
            initial_index = 1;
        end
        %we want to keep the first update persistently, so we shift everything except the first column
        opts.S(:, initial_index:opts.r - 1) = opts.S(:, initial_index + 1:opts.r);
        opts.utility_matrix(:, initial_index:opts.r - 1) = opts.utility_matrix(:, initial_index + 1:opts.r);
        if strcmpi(opts.update, 'bfgs') || strcmpi(opts.update, 'dfp')
            opts.Y(:, initial_index:opts.r - 1) = opts.Y(:, initial_index + 1:opts.r);
        end
        opts.S(:, opts.r) = s_k;

        switch lower(opts.update)
            case 'sr1'
                opts.utility_matrix(:, opts.r) = y_k - A_low_s_k; 
            case {'bfgs', 'dfp'}
                opts.utility_matrix(:, opts.r) = A_low_s_k;
                opts.Y(:, opts.r) = y_k;
        end

    end

end