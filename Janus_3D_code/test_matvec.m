function [ K,KV,Error,KL,ModAll,ModLap] = test_matvec(p,type,Targets,SourceCenters,Targetradii,sourceradii,lambda,V)

np = @(p) 2*(p+1)*p;
sp =@(p) (p+1)^2;



if strcmp('SMat',type)
    params=Setparams(p,SourceCenters,sourceradii,'SL_LMOD_3D','',1,1,3,1,lambda);
    partype='SL_LMOD_3D';
elseif strcmp('DMat',type)
    params=Setparams(p,SourceCenters,sourceradii,'DL_LMOD_3D','',1,1,3,1,lambda);
    partype='DL_LMOD_3D';
elseif strcmp('SpMat',type)
    params=Setparams(p,SourceCenters,sourceradii,'dSL_LMOD_3D','',1,1,3,1,lambda);
    partype='dSL_LMOD_3D';
end
KEparams = Kernel_Eval_parameters(partype,0,1,1,1,1e-8,2,400,1);
KEparams.dim=3;
KEparams.lambda=lambda;
KEparams.W2=params.W2;

KEparams.a=0.5;
params.a=0.5;
targetparams=Setparams(p,Targets,Targetradii,partype,'',1,1,3,1,lambda);
Xtarg=targetparams.X;
Xsource=params.X;
KEparams.nor=params.nor;
if strcmp(type,'SpMat')
    KEparams.nor=targetparams.nor;
K=Kernel_Eval(Xtarg,Xsource,KEparams);

swgt=rand(sp(p),1);


Ntarg=targetparams.nor;


KernelV=K*V;
KV=VSh_Mod_MatVec_RB_trg(V,Xtarg,Ntarg,params);
paramslap=params;
if strcmp(type,'SMat')
    paramslap.flag_pot='SL_L_3D';
elseif strcmp(type,'DMat')
    paramslap.flag_pot='DL_L_3D';
elseif strcmp(type,'SpMat')
    paramslap.flag_pot='dSL_L_3D';
end

KL=VSh_MatVec_RB_trg(V,Xtarg,Ntarg,paramslap);
Error=log10(norm(real(K*V-KV))/norm(real(KV)));


ModAll=VSh_Mod_MatVec_RB2(V,[],params);
ModLap=VSh_MatVec_RB2(V,[],paramslap);


%{
Results(pp,ss)=Error;
    end
end
%}

end

