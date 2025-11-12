function [H, opts] = updateHk(s, y, opts)
switch lower(opts.qn.update)
    case 'bfgs'
        [H, opts] = get_H_BFGS(s, y, opts);
    case 'sr1'
        [H, opts] = get_H_SR1(s, y, opts);
    otherwise
        error([opts.qn.update ' update not implement'])
end
end % updateHk