function M=stokes_smooth_matrix(src,trg)
% src and trg are assumed to be one column in [x;y;z] 
if ~nargin, test_this();return;end

src = reshape(src,[],3);src=src';
trg = reshape(trg,[],3);

ns = size(src,2);
nt = size(trg,1);

M  = zeros(3*nt,3*ns);
P  = repmat(1:3,3,1);Q=P';Q=Q(:);P=P(:);I3=P==Q;
re = reshape((1:3*ns)',[],3)'; re = re(:);

for iT=1:nt
    r(1,:) = trg(iT,1) - src(1,:);
    r(2,:) = trg(iT,2) - src(2,:);
    r(3,:) = trg(iT,3) - src(3,:);
    
    rho = repmat(1./sqrt(dot(r,r,1)),3,1);
    r   = r.*rho;

    Stk = zeros(9,ns);
    Stk(I3,:)=1;
    Stk = Stk+r(P,:).*r(Q,:);
    Stk = Stk.*repmat(rho,3,1);
    Stk = reshape(Stk,3,[]);
    M(iT:nt:end,re) = Stk;    
end
end

function test_this()
src = [0 0 0]';
trg = eye(3);

M = stokes_smooth_matrix(src,trg(:));
R = zeros(3,9);R(:,[1,5,9])=ones(3)+eye(3);R=reshape(R,[],3);
assert(all(all(M==R)));


rng(1984);
ns  = 1000; nt=1000;
src = rand(3*ns,1);trg=2+rand(3*nt,1);den=rand(3*ns,1);

pot = StokesAlltoAll(src, den, trg);
tic;M = stokes_smooth_matrix(src,trg);toc;
disp(max(abs(pot-M*den)));

end