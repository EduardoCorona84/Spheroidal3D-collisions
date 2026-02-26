function [precond, prev] = SpheroidalMS_setprec(ID, n3, precond_type, p, params, prev)
%{
Set preconditioner for the GMRES mobility solve.

Inputs:
    ID      : inverse block(s) for the diagonal preconditioner.
            Either one Nb-by-Nb matrix (shared by all bodies) or a cell array
            with one Nb-by-Nb block per body.
    n3      : number of bodies.
    precond_type : preconditioner type. Only 'bkdiag' is implemented here.

Outputs:
    prec    : function handle that applies the preconditioner (approximate inverse).
    prev    : cached preconditioner data for reuse in later timesteps.
%}
if nargin < 3 || isempty(precond_type)
    precond_type = 'bkdiag';
end

if strcmp(precond_type, 'bkdiag')
    if ~iscell(ID)
        Nb = size(ID,2);
        N = Nb*n3;
        precond = @(V) reshape(ID*reshape(V,Nb,[]),N,[]);
    else
        Nb = size(ID{1},2);
        N = Nb*n3;
        precond = @(V) LOCAL_apply_blockdiag(ID, V, Nb, n3, N);
    end

    prev = ID;
else
    error('Only "bkdiag" preconditioner is implemented for spheroidal_mobility.');
end

end

function Y = LOCAL_apply_blockdiag(blocks, V, Nb, n3, N)
    V = reshape(V, N, []);
    Y = zeros(N, size(V,2));

    for k = 1:n3
        idx = (1:Nb) + Nb*(k-1);
        Y(idx,:) = blocks{k}*V(idx,:);
    end
end
