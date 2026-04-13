function [DMV,D] = RBS_Build_DMV(dtype,Fparams,a,flag_pot)   

if ~strcmp(Fparams.typeMV,'Vsh') || ~(strcmp(Fparams.parslv.prtype,''))

params=Fparams.parbd; p = params.p; 
np=params.np; Sc=params.Sc; eta=Fparams.eta; 
DMV = cell(size(Sc,2),1); D=DMV;

Sc{p}.dblLayerOpStale=1;

%rda = rd;  

for j=1:size(Sc,2)

switch flag_pot    
    case 'TSL_Stk_3D'
    TSSDd=kerneldS_RI(Sc{p},[],1);
     
    if strcmp(dtype,'dense')
        D{j} = a*eye(3*np) + TSSDd; %+Bk'*Ck;       
        DMV{j} = @(x) D{j}*x;
    else
        D{j}=[]; 
    end  

    case 'SL_Stk_3D'
    SSDd=kernelS_RI(Sc{p},[],1);       
     
    if strcmp(dtype,'dense')
        D{j} = SSDd; 
        DMV{j} = @(x) D{j}*x;     
    else
        D{j}=[];  
    end
    
    case 'dSL_L_3D'
    [~, KLDd, ~] = kernelDLap(Sc{p},'SpMat');    
     
    if strcmp(dtype,'dense')
        D{j} = a*eye(np)+eta*KLDd; %+Bk'*Ck;       
        DMV{j} = @(x) D{j}*x;
    else
        D{j}=[]; 
    end      
        
    case 'SL_L_3D'
    [SLDd, ~, ~] = kernelDLap(Sc{p},'SMat');          
     
    if strcmp(dtype,'dense')
        D{j} = SLDd; 
        DMV{j} = @(x) D{j}*x;     
    else
        D{j}=[];  
    end    
end
end

if size(Sc,2)==1
   D = D{1}; 
   DMV = DMV{1}; 
end
else
   D = []; DMV=[];  
end
end