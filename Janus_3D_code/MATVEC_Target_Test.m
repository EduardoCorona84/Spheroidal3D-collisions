clear;clc;
rng('default');
p=16;
C=[0,0,0];
Ctarg=[0,0,0]
lambda=1e-8;
rtarg=[4]';
targparams=Setparams(p,Ctarg,rtarg,'DL_L_3D','',1,1,3,1,lambda);
Xtrg=targparams.X;
r=[1]';


type='SMat';

if strcmp('SMat',type)
    params=Setparams(p,C,r,'SL_L_3D','',1,1,3,1,lambda);
    params.a=1/2;
    paramsmod=params;
    paramsmod.flag_pot='SL_LMOD_3D';
elseif strcmp('DMat',type)
    params=Setparams(p,C,r,'DL_L_3D','',1,1,3,1,lambda);
    params.a=1/2;
    paramsmod=params;
    paramsmod.flag_pot='DL_LMOD_3D';
end

sp = @(p) (p+1)^2;  %(number of harmonic coefficients);
np = @(p) 2*(p+1)*p; %(number of quadrature points);
ii=(1:sp(p))';
nn=floor(sqrt(ii-1)); %indices of n
%generate random seeds


swgt=rand(sp(p),1);

for i=1:sp(p)        %exponentially decaying weights
    swgt(i)=swgt(i)*(1/2)^(floor(sqrt(i-1)));
end
swgt(sp(p-1):end)=0;

swgt(1:end)=0;
swgt(1)=1;
V=shSyn(swgt);
V=repmat(V,size(C,1),1);



%matrix-vector product for laplace
KernelLV=VSh_MatVec_RB_trg(V,Xtrg,[],params);
%kernel for Laplace
KernelL=VSh_MatVec_RB_trg('Mat',Xtrg,[],params);
%matrix-vector product for modified Laplace
KernelV=VSh_Mod_MatVec_RB_trg(V,Xtrg,[],paramsmod);
%kernel for modified Laplace
Kernel=VSh_Mod_MatVec_RB_trg('Mat',Xtrg,[],paramsmod);


%KS=(Sh_Mod_Kernel_Eval_off(swgt,type,0,1,r,[],[],[],lambda));
%KL=(Sh_Mod_Kernel_Eval_off(swgt,type,0,1,r,[],[],[],lambda));

%error(1)=log10(norm(S-KS)/norm(KS));
%error(2)=log10(norm(S-LS)/norm(LS));
%display(error);

Error(1)=log10(norm(real((KernelV-Kernel*V))/norm(real(KernelV))));
Error(2)=log10(norm(real((KernelLV-KernelL*V)))/norm(real(KernelV)));
Error(3)=log10(norm(real(KernelLV-KernelV))/norm(real(KernelV)));


fprintf('\n Rel error ||KernelV-Kernel*V ||/||KernelV|| = %4.4g',Error(1));
fprintf('\n Rel error ||KernelLV - KernelL*V||/||KernelV|| = %4.4g',Error(2)); 
fprintf('\n Rel error ||KernelLV - KernelV||/||KernelV|| = %4.4g',Error(3)); 


fprintf('\n-------------------------------------------------------------\n');









