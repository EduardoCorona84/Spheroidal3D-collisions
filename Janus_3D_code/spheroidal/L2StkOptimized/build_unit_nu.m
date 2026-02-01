function [nu_x, nu_y, nu_z] = build_unit_nu(X_trg, np, ns)
    if isempty(X_trg)
        nu_x = repmat(reshape([1, 0, 0], 1, 3, 1), np, 1, ns);
        nu_y = repmat(reshape([0, 1, 0], 1, 3, 1), np, 1, ns);
        nu_z = repmat(reshape([0, 0, 1], 1, 3, 1), np, 1, ns);
        return;
    end

    nt = size(X_trg, 1);
    nu_x = repmat(reshape([1, 0, 0], 1, 3, 1), nt, 1, ns);
    nu_y = repmat(reshape([0, 1, 0], 1, 3, 1), nt, 1, ns);
    nu_z = repmat(reshape([0, 0, 1], 1, 3, 1), nt, 1, ns);
end