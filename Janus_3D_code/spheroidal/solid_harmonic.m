function F=solid_harmonic(p,u0,u_x)
    F=1;
    if u_x - u0 < 1e-14
        PQ=legendre_otc(p,u_x,0);
        P=PQ{1};
        F=P';
    elseif u_x - u0 > 1e-14
        PQ=legendre_otc(p,u_x,1);
        Q=PQ{2};
        F=Q';
    end
end