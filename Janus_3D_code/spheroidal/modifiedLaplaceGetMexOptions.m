function opts = modifiedLaplaceGetMexOptions(mex_opts)
% Returns defaulted struct or the given options.
    if nargin < 1 || isempty(mex_opts)
        opts = struct('precision_bits', 128, 'output_digits', 17, 'max_memory', 2000);
    else
        error('Invalid MEX options passed in.')
    end
end