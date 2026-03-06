function a = ms_normalization(n, m)
% Meixner-Schafke normalization constant
    a = sqrt(2 / (2*n + 1) * factorial(n + m) / factorial(n - m));
end