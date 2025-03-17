function normal = computeNormal(par, x)
    A = par.R*diag([par.a par.b par.c].^(-2))*par.R.';
    normal = A*(x - par.C);
end