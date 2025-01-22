function Y = RBS_MatVec(V,L,typeMV,params,kerd,a,flag_pot,DMV)

params.flag_pot = flag_pot; params.kerd=kerd; params.a=a; 

if strcmp(typeMV,'Vsh')
    % Build dense matrix
    if isempty(V) && params.dense==1
       V = 'Mat'; 
    end
    
    Y = VSh_MatVec_RB2(V,L,params);
else
    Y = Shg_MatVec(V,DMV,L,params);
end
    
end