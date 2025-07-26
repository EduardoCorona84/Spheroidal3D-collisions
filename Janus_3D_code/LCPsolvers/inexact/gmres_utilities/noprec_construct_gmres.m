function myGmres = noprec_construct_gmres(F, Bk, Lk, TD, SD, Ck, restart, max_iter)

    %This function constructs the function handle for the GMRES solver without preconditioning
    myGmres =@myGmresFunction;

    %Construct some intermediate matrices
    BkF = Bk.'*F;
    FCk = F.'*Ck;

    %Now compute norms used in absolute error criterion
    normF = norm(F', 2);
    normCk = normest(Ck, 2);
    normSD = norm(SD);
    normTDinvs = 1/min(svd(TD));

    function [output, iters, residual] = myGmresFunction(input, tol_abs)

        b = -real(TD*BkF*input) + Lk*BkF*input;
        normB = norm(b);
        tol_rel = tol_abs / (normF * normCk * normSD * normTDinvs * normB);
        [output, ~, residual, iterations] = gmres(TD, b, restart, tol_rel, max_iter);
        residual = residual*normF * normCk * normSD * normTDinvs * normB;
        iters = prod(iterations);
        output = real(FCk*real(SD*(output + BkF*input)));
        
    end
end 