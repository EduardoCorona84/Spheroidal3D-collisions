function A = SHg_Tensor_Kernel_Eval(ind,D,params,type)

nlp = length(params.NP); 

shIdx   = ind(:,1:nlp)-1;
centIdx = ind(:,nlp+1:end)-1;

if size(ind,2)>nlp 
    I = mod(centIdx,2);
    J = (centIdx-I)/2;    
    I = bi2de(I);
    J = bi2de(J);
else
    I = []; J = []; 
end

if nlp>1
    rNP = repmat(params.NP,size(ind,1),1);
    pNP = prod(params.NP); 
    
    IP = mod(shIdx, rNP); 
    JP = (shIdx-IP)./rNP;
    
    pw = repmat(cumprod([1 params.NP(1:end-1)]),size(ind,1),1);
    IP = sum((pw.*IP).').';
    JP = sum((pw.*JP).').';     
else
    pNP = params.NP; 
    IP = mod(shIdx, params.NP); 
    JP = (shIdx-IP)./params.NP;
end

if ~isempty(I)
    Dist = abs(I-J);
    isdiag = ~(Dist>0);    
else
    isdiag = true(size(IP)); 
end
   
I = I*pNP;        
J = J*pNP;  

if ~isempty(I) 
   II = IP+I+1; 
   JJ = JP+J+1;    
else
   II = IP+1; 
   JJ = JP+1;    
end

A = zeros(size(ind,1),1); 

% Diagonal block entries (may change to a sparse matrix format)
if strcmp(type,'dense')
    A(isdiag) = D(IP(isdiag)+1+pNP*JP(isdiag));
elseif sum(isdiag)>0    
    Id = eye(pNP); 
    Id = Id(:,JP(isdiag)+1);      
    Tmp = D(Id);     
    A(isdiag) = Tmp(IP(isdiag)+1+pNP*(0:size(Tmp,2)-1).');    
end

% Off-diagonal block entries
isoffd = ~isdiag; 
if sum(isoffd)>0
params.keval='ind'; 
Aoff = Kernel_Eval_vec(II(isoffd),JJ(isoffd),params);   
A(isoffd) = Aoff;
end