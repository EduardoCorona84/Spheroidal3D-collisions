function [DMV,D] = SpheroidalMS_Build_DMV(dtype,Fparams,a,flag_pot)   
%{
Builds a dense matvec (DMV) for requested operator in flag_pot.

Inputs
dtypes : string
    type of matvec to use (should be 'dense' or empty)
Fparams : struct
    parameters for simulation
a : scalar
    the constant in the operator a*I + O, where O is the layer potential operator
flag_pot : string
    type of layer potential to use; see code below for types:
        'SL_Stk_3D'
        'TSL_Stk_3D'
        'dSL_L_3D'
        'SL_L_3D'

Outputs
DMV : function
    returns the matrix-free operator that represents a*I + O
D : 3*np*n3 x 3*np*n3 matrix
    returns the dense matrix if requested
%}
if ~strcmp(Fparams.typeMV,'SSph') || ~strcmp(Fparams.parslv.prtype,'')
    params=Fparams.parbd; p = params.p; 
    np=Fparams.np; Sc=Fparams.Sc; eta=Fparams.eta; 
    DMV = cell(size(Sc,2),1); D=DMV;
    
    Sc{p}.dblLayerOpStale=1;  
    
    for j=1:size(Sc,2)
        switch flag_pot    
            case 'TSL_Stk_3D'
                TSSDd=kerneldS(Sc{p},[],1);
                if strcmp(dtype,'dense')
                    D{j} = a*eye(3*np) + TSSDd;      
                    DMV{j} = @(x) D{j}*x;
                else
                    D{j}=[]; 
                end  
            case 'SL_Stk_3D'
                SSDd=kernelS(Sc{p},[],1);       
                if strcmp(dtype,'dense')
                    D{j} = SSDd; 
                    DMV{j} = @(x) D{j}*x;     
                else
                    D{j}=[];  
                end
            case 'dSL_L_3D'
                [~, KLDd, ~] = kernelDLap(Sc{p},'SpMat');    
                if strcmp(dtype,'dense')
                    D{j} = a*eye(np)+eta*KLDd;   
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