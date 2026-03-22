function opts = modified_laplace_get_mex_options(mex_opts)
    if nargin < 1 || isempty(mex_opts)
        mex_opts = struct();
    end

    opts = mex_opts;
    if ~isfield(opts, 'precision_bits') || isempty(opts.precision_bits)
        opts.precision_bits = 128;
    end
    if ~isfield(opts, 'output_digits') || isempty(opts.output_digits)
        opts.output_digits = 17;
    end
    if ~isfield(opts, 'max_memory') || isempty(opts.max_memory)
        opts.max_memory = 2000;
    end
end
