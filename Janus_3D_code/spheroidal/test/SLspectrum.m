function [lambda_int,lambda_surf,lambda_ext]=SLspectrum(p,u0,a)
    sp=(p+1)^2;
    ii = (1:sp)'; nn = floor(sqrt(ii-1)); mm=ii-nn.^2-nn-1;
    
    cnm=a.*factorial(nn-mm)./factorial(nn+mm).*(-1).^mm .*sqrt(u0.^2-1);
    L=legendre_otc(p,u0,1);
    P=L{1}; Q=L{2};

    lambda_int=cnm.*Q;
    lambda_surf=cnm.*P.*Q;
    lambda_ext=cnm.*P;   
end