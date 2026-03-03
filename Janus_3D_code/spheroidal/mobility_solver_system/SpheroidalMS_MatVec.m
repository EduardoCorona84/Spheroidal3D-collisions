function Y = SpheroidalMS_MatVec(V,L,typeMV,Fparams,kerd,a,flag_pot)
%{
Intermediatery function to handle processing to different mavecs.

Inputs
    V : density (and so, [] passed in means build the matrix)
    L : nullspace correction/completion term
    typeMV : type of matvec; see spheroidal_mobility.m
    Fparams : mobility solver params
    kerd : kernel dimension (should probably always be 3)
    a : the scalar in the BIE a*I + O, where O is the operator dictated by flag_pot
    flag_pot : the layer potential to compute
    DMV : 
Outputs
    Y : the BIE a*I + A, where A is the operator dictated by flag_pot
%}
 
Fparams.flag_pot = flag_pot;
Fparams.kerd=kerd;
Fparams.a=a; 

if isempty(V) && Fparams.dense==1
    V = 'Mat'; 
end

if strcmp(typeMV,'SSph')
    Y = SSph_MatVec(V,L,Fparams);
elseif strcmp(typeMV, 'vectorSSph')
    error('Not implemented.');
else
    error("This mode of matvec is not implemented.");
end
end
