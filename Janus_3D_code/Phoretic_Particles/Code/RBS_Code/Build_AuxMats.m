function [C,B,D,L] = Build_AuxMats(Wg,Xg,Xc,np,n3)

N = 3*np*n3; 
IC = []; JC = []; VC = []; 
VD = []; VB = [];  

if nargout==4
IL = []; JL = []; VL = []; 
end

for i=1:n3
    
idx = (1:np)+np*(i-1); 
X = Xg(idx,:);
W = Wg(idx); 

oW = ones(size(W));
sumW = sum(W); 

W2 = zeros(1,3*np); 
W2(1:3:3*np)=W; W2(2:3:3*np)=W; W2(3:3:3*np)=W; 
    
if ~isempty(Xc)
% Center X
X = X - repmat(Xc(i,:),np,1); 
end

indx =(1:3:3*np)+3*np*(i-1); 
indy =(2:3:3*np)+3*np*(i-1); 
indz =(3:3:3*np)+3*np*(i-1);

% C*sigma = [int{sigma} ; int{X \times sigma}]    
IC = [IC; reshape(repmat(6*(i-1)+(1:6),np,1),[],1) ; reshape(repmat(6*(i-1)+(4:6),np,1),[],1)];
JC = [JC; indx'; indy'; indz'; indy'; indz'; indx'; indz'; indx'; indy']; 
VC = [VC; W; W; W ; -W.*X(:,3); -W.*X(:,1); -W.*X(:,2); W.*X(:,2);  W.*X(:,3); W.*X(:,1)];
% First three vectors are just integrals of f for each coordinate
VD = [VD; oW; oW; oW ; -X(:,3); -X(:,1); -X(:,2); X(:,2);  X(:,3); X(:,1)];

tau1 = sum(W.*X(:,3).^2)+sum(W.*X(:,2).^2); 
tau2 = sum(W.*X(:,1).^2)+sum(W.*X(:,3).^2); 
tau3 = sum(W.*X(:,2).^2)+sum(W.*X(:,3).^2); 

VB = [VB; oW./sumW; oW./sumW; oW./sumW ; ...
    -X(:,3)./tau1; -X(:,1)./tau2; -X(:,2)./tau3; ...
     X(:,2)./tau1;  X(:,3)./tau2; X(:,1)./tau3];

end

C = sparse(IC,JC,VC,6*n3,N); 
B = sparse(IC,JC,VB,6*n3,N);
D = sparse(IC,JC,VD,6*n3,N);

if nargout==4
   for i=1:n3
      % We divide by the integral of ((X-X_c) (x) 1)^2 over Sc. 
      idxM = (1:3*np)+3*np*(i-1); 
      [JJ,II] = meshgrid(idxM,idxM); 
      IL = [IL ; II(:)]; 
      JL = [JL ; JJ(:)]; 
      Lb = B(6*(i-1)+(1:6),idxM)'*C(6*(i-1)+(1:6),idxM);   
      VL = [VL ; Lb(:)]; 
   end
   L = sparse(IL,JL,VL,N,N);  
end

end