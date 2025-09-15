function opts = get_correction_SR1(k, fg_low, s_k, y_k, opts, info)

    %TO DO: implement update storage and compact representation
    %Three cases:
    %1. Dense Full Memory (with and without skips)
    %2. Dense Limited Memory (with and without skips)
    %3. Compact Representation (no skips for now, want to add)

    switch lower(opts.memory)
        case 'dense full'
            % Implement dense full memory update, we just store a full dense matrix that's accumulating all the update. Have an option for skips
            [~, ~, A_low_s_k] = fg_low(s_k);
            quantity = y_k - (A_low_s_k + opts.update_matrix*s_k);
            if opts.skips == true
                %check skip condition and skip if so
                if abs(s_k'*quantity) < 1e-8*norm(s_k)*norm(quantity)
                    %return and do not update opts.outer.update_matrix
                    return
                end
            end
            %If no skip, compute quantity of the needed matrix
            opts.update_matrix = opts.update_matrix + ((quantity) * (quantity)' / (s_k' * quantity));

        case 'dense limited'
            % Implement dense limited memory update, we store r amount of y_k and s_k and then use these for the update. The update will be a full dense matrix, but this lets us use limited memory. Have an option for skips.

        case 'compact'
            % Implement compact representation update. Not sure how to handle skips, if at all, yet.

        otherwise
            error('Unknown memory option');
    end
end


