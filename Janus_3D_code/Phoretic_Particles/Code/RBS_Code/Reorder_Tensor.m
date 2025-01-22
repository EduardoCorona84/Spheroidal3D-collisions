function A = Reorder_Tensor(A,p,type,J)

NP1 = factor(p+1);     L1 = length(NP1); 
NP2 = factor(p); NP2 = NP2(end:-1:1); NP2 = [NP2 2]; L2 = length(NP2);

if nargin<4   
I2 = L2:2:L1+L2; 
I1 = 1:L1+L2; I1(I2)=[];
J = zeros(L1+L2,1); 
J(I1) = 1:L1; 
J(I2) = L1+1:L1+L2;
end

if isempty(A)
   NP = [NP1 NP2];
   A = NP(J);    
else
if strcmp(type,'vec')
NP = [NP1 NP2 size(A,2)];
sA = size(A); 
A = reshape(permute(reshape(A,NP),[J ; (L1+L2+1)]),sA);    
else
NP = [NP1 NP2 NP1 NP2];
J = [J J+L1+L2]; 
sA = size(A); 
A = reshape(permute(reshape(A,NP),J),sA);
end
end
   
end