function psi = constructPsiGlobal(theta)
    P = [0 -theta(4) theta(3);
    theta(4) 0 -theta(2);
    -theta(3) theta(2) 0];

    psi = (1/2).*[-theta(2:4).' ; 
             theta(1).*eye(3) - P];
end