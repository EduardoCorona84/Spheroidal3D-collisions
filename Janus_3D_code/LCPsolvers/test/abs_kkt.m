function e = abs_kkt(x, A, b)
phi = min(x,A*x + b);
e = dot(phi, phi);
end % abs_kkt