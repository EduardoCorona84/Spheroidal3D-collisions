function opts = get_correction_DFP(k, fg_low, s_k, y_k, opts, info)

    %TO DO: implement update storage and compact representation
    %Three cases:
    %1. Dense Full Memory (with and without skips)
    %2. Dense Limited Memory (with and without skips)
    %3. Compact Representation (no skips for now, want to add)

    switch lower(opts.memory)
        case 'dense full'
            % Implement dense full memory update, we just store a full dense matrix that's accumulating all the update. Have an option for skips
            [~, ~, A_low_s_k] = fg_low(s_k);
            

        case 'dense limited'
            % Implement dense limited memory update, we store r amount of y_k and s_k and then use these for the update. The update will be a full dense matrix, but this lets us use limited memory. Have an option for skips.

            %no skip implementation for now, very simple


        case 'compact'

            %no skip implementation for now, very simple. 
            %need to store S and Y - B_0*S matrices, we will just treat y_k - B_0*s_k as y_k for storage purposes
            [~, ~, A_low_s_k] = fg_low(s_k);
            opts = update_correction_memory(k, s_k, y_k, A_low_s_k, opts);
            %using the multisecant form, this may have instabiliy problems, see inner_solver for gradient formulation to see the evaluation.

        
        otherwise
            error('Unknown memory option');
    end
end
