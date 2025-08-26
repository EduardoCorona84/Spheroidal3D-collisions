function Y = RBS_MatVec(V,L,typeMV,params,kerd,a,flag_pot,DMV)

params.flag_pot = flag_pot; params.kerd=kerd; params.a=a; 

if strcmpi(typeMV,'Vsh')
    % Build dense matrix
    if isempty(V) && params.dense==1
       V = 'Mat'; 
    end
 if(isfield(params,'lambda'))   
    Y=VSh_Mod_MatVec_RB2(V,L,params);
 else
     Y = VSh_MatVec_RB2(V,L,params);
 end
else
    Y = Shg_MatVec(V,DMV,L,params);
end
    
end