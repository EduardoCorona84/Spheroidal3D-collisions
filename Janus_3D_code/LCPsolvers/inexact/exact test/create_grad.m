function gradfunc = create_grad(A,b);

gradfunc = @fcnGrad;

function [fk, gradk] = fcnGrad(x)
    fk = 0.5 * dot(x, A * x) + dot(b, x);
    gradk = A * x + b;
end

end