% misc extra parameters
tol=1e-4; 
mdist=3; 
denseMV=false; 
denseforce=1;
gamma=1; 
parslv = struct('solver','gmres','tol',tol,'maxit',200,'rst',4,'prtype','bkdiag','prec',[],'prLCP',false);  
ix = 2
F = A_list{2}.F;
Ck = A_list{2}.Ck;
SD = A_list{2}.SD;
TD = A_list{2}.TD;
Bk = A_list{2}.Bk;
Lk = A_list{2}.Lk;
Bf = @(x) (Bk.')*(F*x);

A = @(x) real(F.'*(Ck*Lapp(SD,Lslv(TD,-Lapp(TD,Bf(x))+Lk*Bf(x),parslv)+Bf(x))));