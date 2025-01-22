function [SMat, SpMat, DMat] = kernelModifiedLap(S,type,lambda)
  persistent p Ywt A Wsph    %S is sphere shape

  if ( nargin > 0 && ~isempty(S) )
    printMsg('  * Generating the direct double-layer matrix for p=%d.\n',S.p);
   
    np = length(S.cart.x);
    if (isempty(p) || S.p~=p)
      [Ywt A Wsph] = getPesistents(S.p);
      p = S.p;
    end
  
    W0 = S.geoProp.W./Wsph;
    DMat = zeros(np, np); SMat = DMat; SpMat = DMat;
    ind_ref = reshape((1:np),p+1,2*p);

    for k = 1:2*p
      ind = circshift(ind_ref,[0 k-1]);
      ind = ind(:);
      
      for j = 1:p+1
        %- Rotating the pole to (j,k)
        R = A{j}(ind,ind);
        xx = (R*S.cart.x)'; nx = (R*S.geoProp.nor.x)';
        yy = (R*S.cart.y)'; ny = (R*S.geoProp.nor.y)';
        zz = (R*S.cart.z)'; nz = (R*S.geoProp.nor.z)';
        
        W = ((R*W0).*Wsph)';
        ind_pole = (k-1)*(p+1) + j;

        %- The distance vector
        rx = S.cart.x(ind_pole) - xx;
        ry = S.cart.y(ind_pole) - yy;
        rz = S.cart.z(ind_pole) - zz;
        
        rho=sqrt(rx.^2 + ry.^2 + rz.^2);
        inv_rho=1./rho;
        %- Normalizing
        rx = rx.*inv_rho;
        ry = ry.*inv_rho;
        rz = rz.*inv_rho;

        %- Laplace kernel and the scalar part
        g = (nx.*rx + ny.*ry + nz.*rz);
        %g = g.*Ywt.*W.*inv_rho.*inv_rho;

       
       
            switch type
                case 'SMat'
                 %- Single Layer Matrix Modified Laplace
                SMat(ind_pole,:)=(Ywt.*W.*inv_rho.*exp(-lambda*rho))*R;          
                case 'SpMat'
                %- Derivative of SLP
                SpMat(ind_pole,:) = -((S.geoProp.nor.x(ind_pole)*rx + S.geoProp.nor.y(ind_pole)*ry + S.geoProp.nor.z(ind_pole).*rz).*Ywt.*W.*inv_rho.*inv_rho)*R;    
                case 'DMat'
                %- Double layer Matrix Modified Laplace
                
                g=g.*Ywt.*W;
                Expon=exp(-lambda*rho);
                Paren=(inv_rho.*inv_rho +lambda*inv_rho);
                DMat(ind_pole,:)=(g.*Expon.*Paren)*R;
                
               % DMat(ind_pole,:)=((exp(-lambda*rho).*inv_rho).*g.*Ywt.*W.*(inv_rho.*inv_rho+lambda*inv_rho))*R;
               
               
            end
        end
        
%         f_Temp = -((S.geoProp.nor.x(ind_pole)*rx + S.geoProp.nor.y(ind_pole)*ry + S.geoProp.nor.z(ind_pole).*rz).*W.*inv_rho.*inv_rho);
%         surf(reshape((g./Ywt)',p+1,2*p))
%         pause
      end
    end
    
  

function [Ywt, A, Wsph] = getPesistents(p)

  printMsg('  * Generating the direct double-layer integration data for p=%d.\n',p);
  Ywt = (1/4/pi)*SingularWeights(p);
  u = gl_grid(p);
  Wsph = sin(u);
  np = 2*p*(p+1);

  A = cell(p+1,1);
  for idx = 1:p+1
    fname = ['RotMat' num2str(idx) '-'];
    A{idx} = readData(fname, p, [np np]);
    
    if(isempty(A{idx}))
      printMsg('  * Generating the %dth direct rotation matrices for p=%d.\n',idx,p-1);
      A{idx} = movePole(eye(np),u(idx),0);
      writeData(fname, p, A{idx});
    end
  end
  printMsg('  * Direct rotation matrices for p=%d were generated/read from file.\n',p);
end

end

