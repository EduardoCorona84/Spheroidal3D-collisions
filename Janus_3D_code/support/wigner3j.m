function w = wigner3j(l1, l2, l3, m1, m2, m3)
%{
Function to calculate Wigner-3j via the DLMF formula.   
%}
[l1, l2, l3, m1, m2, m3] = local_broadcast(double(l1), double(l2), double(l3), double(m1), double(m2), double(m3));
w = zeros(size(l1));

l1 = round(l1); l2 = round(l2); l3 = round(l3);
m1 = round(m1); m2 = round(m2); m3 = round(m3);

% Selection rule
valid = (m1 + m2 + m3 == 0);
valid = valid & (abs(m1) <= l1) & (abs(m2) <= l2) & (abs(m3) <= l3);
valid = valid & (l3 >= abs(l1 - l2)) & (l3 <= l1 + l2);
valid = valid & (mod(l1 + l2 + l3, 2) == 0);

if ~any(valid, 'all')
    return;
end

nmax = max(l1(valid) + l2(valid) + l3(valid)) + 1;
logfact = local_logfact_cache(nmax);

idx = find(valid);
for k = 1:numel(idx)
    i = idx(k);
    w(i) = local_wigner3j_scalar(l1(i), l2(i), l3(i), m1(i), m2(i), m3(i), logfact);
end
end


function w = local_wigner3j_scalar(l1, l2, l3, m1, m2, m3, logfact)
    kmin = max([0, (l2 - l3 - m1), (l1 + m2 - l3)]);
    kmax = min([(l1 + l2 - l3), (l1 - m1), (l2 + m2)]);
    if kmin > kmax
        w = 0;
        return;
    end

    % Delta
    a = l1 + l2 - l3;
    b = l1 - l2 + l3;
    c = -l1 + l2 + l3;
    d = l1 + l2 + l3 + 1;
    logDelta = 0.5 * (logfact(a + 1) + logfact(b + 1) + logfact(c + 1) - logfact(d + 1));

    % sqrt((l1+m1)!(l1-m1)!...(l3-m3)!).
    logSqrtFac = 0.5 * ( ...
        logfact(l1 + m1 + 1) + logfact(l1 - m1 + 1) + ...
        logfact(l2 + m2 + 1) + logfact(l2 - m2 + 1) + ...
        logfact(l3 + m3 + 1) + logfact(l3 - m3 + 1));

    % (-1)^(l1-l2-m3).
    phase = 1 - 2 * mod(l1 - l2 - m3, 2);

    k = kmin:kmax;
    logDen = logfact(k + 1) + ...
             logfact(a - k + 1) + ...
             logfact(l1 - m1 - k + 1) + ...
             logfact(l2 + m2 - k + 1) + ...
             logfact(l3 - l2 + m1 + k + 1) + ...
             logfact(l3 - l1 - m2 + k + 1);
    alt = 1 - 2 * mod(k, 2);  % (-1)^k
    s = sum(alt .* exp(-logDen));

    w = phase * exp(logDelta + logSqrtFac) * s;
end


function [a, b, c, d, e, f] = local_broadcast(a, b, c, d, e, f)
    sz = [size(a); size(b); size(c); size(d); size(e); size(f)];
    is_scalar = prod(sz, 2) == 1;
    non_scalar = find(~is_scalar, 1, 'first');
    if isempty(non_scalar)
        return;
    end
    target = sz(non_scalar, :);
    if is_scalar(1), a = repmat(a, target); end
    if is_scalar(2), b = repmat(b, target); end
    if is_scalar(3), c = repmat(c, target); end
    if is_scalar(4), d = repmat(d, target); end
    if is_scalar(5), e = repmat(e, target); end
    if is_scalar(6), f = repmat(f, target); end
end


function logfact = local_logfact_cache(nmax)
    persistent LF NMAX
    if isempty(LF)
        NMAX = -1;
    end
    if nmax > NMAX
        NMAX = nmax;
        LF = gammaln((0:NMAX) + 1); % LF(n+1)=log(n!)
    end
    logfact = LF;
end
