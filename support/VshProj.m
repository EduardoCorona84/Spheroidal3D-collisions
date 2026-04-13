function PF = VshProj(F,type)

NF = size(F,1); 
PF = zeros(size(F)); 
J = [1:3:NF 2:3:NF 3:3:NF]; 
PF(J,:) = VshSyn(VshAna(F(J,:),type),type); 

end