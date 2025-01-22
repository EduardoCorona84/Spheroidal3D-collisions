p=8;
np = 2*p*(p+1);

pars = SpheroidalParameters;
pars.u0=1.2;
pars.a=1;
pars.sigma = eye(np);

X = prolate_spheroid_shape(p,pars.u0,pars.a);
Xeval = [X+1.2;X-1.2];

pars.plot(Xeval)

DMoff = spheroidalDL(pars,Xeval);


% Compare with known Ynm
sigma=Yp(2,1,p);
DMsigma_off = DMoff*sigma;

pars.sigma=sigma;
Dsigma_off = spheroidalDL(pars,Xeval);

err=abs(DMsigma_off-Dsigma_off);
figure;
semilogy(1:length(err),err);

% ---------------------------------------------
p=8;

params=SpheroidalParameters;
params.u0 = [1.1 1.2];
params.a=1./params.u0;
params.centers=[0 0 0; 1 1 1];
params.thetas = [0 pi/3];
params.phis=[0 pi/5];
params.matvec_eta = 10;

M=spheroidalMatVecKernel(params,'DL',p);

sigma=permute([Yp(2,1,p) Yp(2,2,p)],[1 3 2]);
params.sigma = sigma;

DMsigma=M*sigma(:);
Dsigma= spheroidalMatVec(params,'DL');
Dsigma=Dsigma(:);

err=abs(Dsigma-DMsigma);
max(abs(err))


