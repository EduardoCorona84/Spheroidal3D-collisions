function [lambda_int,lambda_surf,lambda_ext]=DLspectrum(p,u0)
    sp=(p+1)^2;
    ii = (1:sp)'; nn = floor(sqrt(ii-1)); mm=ii-nn.^2-nn-1;
    anm=factorial(nn-mm)./factorial(nn+mm).*(-1).^mm .*(u0.^2-1);
    L=legendre_otc(p,u0,1,1,1);
    P=L{1}; Q=L{2}; dP=L{3}; dQ=L{4};

    lambda_int=anm.*dQ;
    lambda_surf=anm./2.*(P.*dQ+dP.*Q);
    lambda_ext=anm.*dP;   
end