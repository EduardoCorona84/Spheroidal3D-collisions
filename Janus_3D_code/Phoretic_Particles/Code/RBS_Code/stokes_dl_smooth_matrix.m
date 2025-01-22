function M=stokes_dl_smooth_matrix(src,nor,trg)
% src and trg are assumed to be one column in [x;y;z] 
if ~nargin, test_this();return;end

src = reshape(src,[],3);src=src';
nor = reshape(nor,[],3);nor=nor';

trg = reshape(trg,[],3);

ns = size(src,2);
nt = size(trg,1);

M  = zeros(3*nt,3*ns);
P  = repmat(1:3,3,1);Q=P';Q=Q(:);P=P(:);
re = reshape((1:3*ns)',[],3)'; re = re(:);

for iT=1:nt
    r(1,:) = src(1,:) - trg(iT,1);
    r(2,:) = src(2,:) - trg(iT,2);
    r(3,:) = src(3,:) - trg(iT,3);
    
    coeff = sqrt(dot(r,r,1)).^-5;
    coeff = coeff.*dot(r,nor,1);

    Stk = r(P,:).*r(Q,:);
    Stk = Stk.*repmat(coeff,9,1);
    Stk = reshape(Stk,3,[]);
    M(iT:nt:end,re) = Stk;    
end
end

function test_this()
src = [0 0 0]'; nor = [1 1 1]/sqrt(3);
trg = eye(3);

M = stokes_dl_smooth_matrix(src,nor,trg(:));
R = zeros(3,9);R(:,[1,5,9])=-1/sqrt(3)*eye(3);R=reshape(R,[],3);
assert(all(all(M==R)));

rng(1984);
ns  = 1000; nt=1000;
src = rand(3*ns,1);
nor = rand(3*ns,1);
den = rand(3*ns,1);
trg = 2+rand(3*nt,1);
 
pot = DoubleAlltoAll(src, den, nor, trg);
tic;M = stokes_dl_smooth_matrix(src,nor,trg);toc;
disp(max(abs(pot-M*den)));

end