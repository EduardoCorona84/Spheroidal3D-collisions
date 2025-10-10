function Y = SpheroidalMS_MatVec(V,L,typeMV,Fparams,kerd,a,flag_pot,DMV)
%{
Intermediatery function to handle processing to different mavecs.

Inputs
    V : density (and so, [] passed in means build the matrix)
    L : nullspace correction/completion term
    typeMV : type of matvec; see spheroidal_mobility.m
    Fparams :
    kerd : kernel dimension (should probably always be 3)
    a : the scalar in the BIE a*I + O, where O is the operator dictated by flag_pot
    flag_pot : the layer potential to compute
    DMV : 
Outputs
    Y : the BIE a*I + O, where O is the operator dictated by flag_pot
%}
 
Fparams.flag_pot = flag_pot;
Fparams.kerd=kerd;
Fparams.a=a; 

if strcmp(typeMV,'SSph')
    % Build dense matrix
    if isempty(V) && Fparams.dense==1
        V = 'Mat'; 
    end

    if(isfield(Fparams,'lambda')) % Modified Laplace
        error("This problem type is not implemented.");
        Y = VSh_Mod_MatVec_RB2(V,L,Fparams);
    else
        Y = SSph_MatVec(V,L,Fparams);
    end
else % Rotation-based singular quadrature matvec
    error("This mode of matvec is not implemented.");
    Y = Shg_MatVec(V,DMV,L,Fparams);
end
        
end