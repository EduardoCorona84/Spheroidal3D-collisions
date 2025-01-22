function [Df,mvH] = LkernelD(f, S, mvH, directIn,type)   

  persistent matVecHandle lastSurf
  
  if(nargin==0), testKernelD(); return; end
   
  if nargin>2
  if ~isempty(mvH) 
      matVecHandle=mvH; 
  end
  end
  
  if(nargin<5), type='rotfree'; end
   
  pDirect = 100; %empirical    
  direct = max(S.p) < pDirect;
 
  if(nargin>3), direct = directIn; end
  
  buildmatVec = isempty(matVecHandle) || ( any([S.dblLayerOpStale]) || any(S ~= lastSurf)); 
  
  if buildmatVec
  if direct 
    lastSurf = S;
    for ii=1:length(S)
        DMat{ii} = kernelDMatrix(S(ii),type);
        S(ii).dblLayerOpStale = false;
    end
    matVecHandle = @(den) multivec(DMat,den);
  else
    matVecHandle = @(den) LkernelDMatFree(den, S );   
  end
  end   
  
  if nargout>1 && buildmatVec
     mvH = matVecHandle;
  end
  
  Df = matVecHandle(f);
end
    
function v=multivec(A,b)
l = length(A);
v{1:l} = zeros(size(b{1}));

for ii=1:l
    v{ii} = A{ii}*b{ii};
end  
end

% We take kernelDMatrix from Stokes and modify it. Since it computes the
% Laplace Double Layer as part of the Stokes kernel, we use that. 
function DMat = kernelDMatrix(S,type)
  persistent p Ywt YwtIm A Wsph vRot uRot Rotmv

  if ( nargin > 0 && ~isempty(S) )
    %printMsg('* Generating the direct double-layer matrix for p=%d.\n',S.p);
   
    np = length(S.cart.x);
    %if (isempty(p) || S.p~=p)
      [Ywt,A,Wsph] = getPersistents(S.p,type);      
      if strcmp(type,'rotfree')
        Ywt = Ywt.';
        YwtIm = flipud(reshape(Ywt,S.p+1,[])); YwtIm = YwtIm(:);    
        [~,vRot,uRot]=RotMat_apply(S.p,vRot,uRot,[],[]);  
        Rotmv = @(shc,idx) RotMat_apply(S.p,vRot,uRot,shc,idx);   
      else
        vRot=[]; uRot=[]; Rotmv=[];    
      end  
      p = S.p;
    %end  
  
    % DMat is now np x np. 
    DMat = deal(zeros(np));
    
if strcmp(type,'rotfree')
    
    [~, gwt]=g_grid(p+1); 
    wt = pi/p*repmat(gwt',1,2*p);   
    wt = repmat(wt(:),1,2*p);  
    
    dot3 = @(X,Y) X(1:end/3,:).*Y(1:end/3,:) + ...
         X(end/3+1:2*end/3,:).*Y(end/3+1:2*end/3,:) + ...
         X(2*end/3+1:end,:).*Y(2*end/3+1:end,:);
    
    X = reshape(S.cart.to_array,[],3);  
    N = reshape(S.geoProp.nor.to_array,[],3);  
    W = S.geoProp.W./Wsph;       
    
    shc = shAna([X W N]);
    
    %X,W,f index in the rotated matrix
    nCol = 7; 
    ii = reshape(1:nCol*2*p,nCol,[]);
    xInd = ii(1:3,:); xInd = xInd(:);
    wInd = ii(4,:); wInd = wInd(:);
    nInd = ii(5:7,:); nInd = nInd(:);
    
    ii = reshape(1:np,p+1,[]);
    northHem = ii(1:floor(p/2)+1  ,:); 
    southHem = ii(floor(p/2)+2:end,:); 
    southHem = flipud(circshift(southHem,[0 -p]));
 
    for ii=1:floor(p/2)+1
      %Rotating to the desired lattitude
      rot = Rotmv(repmat(shc,1,2*p),ii);   
      %XP is the pole positions (north pole)
      XP = X(northHem(ii,:),:);
      XP = XP.'; XP = XP(:).';
      %The distance in the physical space
      r = rot(:,xInd);
      r = r - repmat(XP,np,1);
      r = reshape(r,3*np,2*p);
      %Rotated surface element times the integration weights
      Wr = rot(:,wInd).*repmat(Ywt.*Wsph,1,2*p);  
      %Rotated normal
      nr = rot(:,nInd);
      nr = reshape(nr,3*np,2*p);  
      %The distance and stokes kernel 
      rho = sqrt(dot3(r,r));
      rdotn = dot3(r,nr); 
      g = Wr.*rdotn./rho.^3; 
      g = RotTrans_apply(p,vRot,uRot,shSynT(g),ii);   
      
      %for k=1:np
      %   ek = zeros(np,1); ek(k)=1; 
      %   ekr = Rotmv(repmat(shAna(ek),1,2*p),ii); 
      %   DMat(northHem(ii,:),k) = sum(g.*ekr).'; 
      %end
      %dumping the results   
      DMat(northHem(ii,:),:) = g.';         
      
      %Treating the south pole at the same time
      if(ii~=p/2+1)
        %XP is the pole positions for the south pole
        XP = X(southHem(ii,:),:);
        XP = XP'; XP = XP(:)';
        %The distance in the physical space
        r = rot(:,xInd);
        r = r - repmat(XP,np,1);
        r = reshape(r,3*np,2*p);
        %Rotated surface element times the integration weights
        Wr = rot(:,wInd).*repmat(YwtIm.*Wsph,1,2*p);      
        %The distance and stokes kernel
        rho = sqrt(dot3(r,r));
        rdotn = dot3(r,nr); 
        g = Wr.*rdotn./rho.^3;
        g = RotTrans_apply(p,vRot,uRot,shSynT(g),ii);        
       
      %for k=1:np
      %   ek = zeros(np,1); ek(k)=1; 
      %   ekr = Rotmv(repmat(shAna(ek),1,2*p),ii);
      %   DMat(southHem(ii,:),k) = sum(g.*ekr).'; 
      %end  
        %dumping the results        
        DMat(southHem(ii,:),:) = g.';                  
      end
    end
    
else
    W0 = S.geoProp.W./Wsph;  
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
        rx = xx - S.cart.x(ind_pole);
        ry = yy - S.cart.y(ind_pole);
        rz = zz - S.cart.z(ind_pole);
        inv_rho = 1./sqrt(rx.^2 + ry.^2 + rz.^2);

        %- Normalizing
        rx = rx.*inv_rho;
        ry = ry.*inv_rho;
        rz = rz.*inv_rho;

        %- Laplace kernel (Ywt*Area element*(n'r)/rho^3
        g = (nx.*rx + ny.*ry + nz.*rz);
        g = g.*Ywt.*W.*inv_rho.*inv_rho;
        
        %
        DMat(ind_pole,:) = g*R;  
      end
    end 
    %}
end
  end
end

function [Ywt,A,Wsph] = getPersistents(p,type)

  %printMsg('* Generating the direct double-layer integration data for p=%d.\n',p);
  Ywt = (1/4/pi)*SingularWeights(p);    
  u = gl_grid(p);
  Wsph = sin(u);
  np = 2*p*(p+1);

  if ~strcmp(type,'rotfree')
  A = cell(p+1,1);
  for idx = 1:p+1
    fname = ['RotMat' num2str(idx) '-'];
    A{idx} = readData(fname, p, [np np]);
    
    if(isempty(A{idx}))
      %printMsg('* Generating the direct rotation matrices for p=%d, idx=%d.\n',p,idx-1);
      A{idx} = movePole(eye(np),u(idx),0);
      writeData(fname, p, A{idx});
    end
  end
  %printMsg('* Direct rotation matrices for p=%d were generated/read from file.\n',p);
  else
     A = []; 
     %printMsg('* Rotation matrices for p=%d are matrix free.\n',p);
  end
end

function [rotx,vRot,uRot] = RotMat_apply(p,vRot,uRot,shc,ui) 

nCol=size(shc,2)/(2*p);      
[~, v] = gl_grid(p);
v = reshape(v,p+1,[]);
v = repmat(v(1,:),nCol,1); v = v(:)';
%Rotation matrix in the longitude direction
if isempty(vRot)
vRot = 1i*repmat(-p:p,p+1,1);
vRot = exp(shrinkShVec(vRot(:))*v);
end
%Rotation matrix in the lattitude direction
if isempty(uRot)
uRot = getMat(p);
end

if ~isempty(shc)
    %Making copies for each v=const and aligning the pole with all of these.
    shc= vRot.*shc;  
    
    %Rotating to the desired lattitude
    shcRot = zeros((p+1)^2,size(shc,2)); 
    for jj = 0:p
        ind = jj^2+1:(jj+1)^2;
        %size(uRot{ui}{jj+1})
        %size(shc(ind,:))
        shcRot(ind,:) = uRot{ui}{jj+1}*shc(ind,:);
    end
    
    shcRot = conj(vRot).*shcRot;   
    
    rotx = shSyn(shcRot, false);
else
    rotx=[]; 
end

end  

function rotx = RotTrans_apply(p,vRot,uRot,shc,ui) 
    %Making copies for each v=const and aligning the pole with all of these.
    [~, v] = gl_grid(p);
    v = reshape(v,p+1,[]);
    v = v(1,:); v = v(:)'; 
    
    %Rotation matrix in the longitude direction
    vRotT = -1i*repmat(-p:p,p+1,1);  
    vRotT = exp(shrinkShVec(vRotT(:))*v);  
    
    shc = conj(vRotT).*shc;   
    
    %Rotating to the desired lattitude
    shcRot = zeros((p+1)^2,size(shc,2)); 
    for jj = 0:p
        ind = jj^2+1:(jj+1)^2;
        shcRot(ind,:) = uRot{ui}{jj+1}.'*shc(ind,:);              
    end
    
    shcRot= vRotT.*shcRot;   
    
    rotx = shAnaT(shcRot, false);  
end  

function Rall = getMat(p)

  nmat = floor(p/2)+1;
  Rall = cell(nmat,1);
  matSize =@(n) n*(4*n^2-1)/3;
  
  pstr = ceil(abs(log10(nmat)));
  pstr = ['%0' num2str(pstr) 'd'];
  for nt=1:nmat 
    R = [];
    R = cell(p+1,1);
    fileName = ['rotMatSh_' num2str(p) '_'];
    count = matSize(p+1); 

    RM = readData(fileName, num2str(nt, pstr), 2*count);
    if(~isempty(RM))
      rReal = RM(1:count);
      rImag = RM(count+1:end);
      for n=0:p
        ind = matSize(n);
        dim = (2*n+1);
        R{n+1} = reshape(complex(rReal(ind+1:ind+dim^2), ...
                                 rImag(ind+1:ind+dim^2)),dim,dim);
      end
    else
      u = gl_grid(p);
      theta = u(nt);
      %printMsg('* Generating the fast rotation data for p=%d and target=%d.\n',p,nt);

      shcIn = zeros((p+1)*(2*p+1),2*p+1);
      for m=1:2*p+1
        ind = (m-1)*(p+1)+1;
        shcIn(ind:ind+p,m) = 1;
      end
      
      fIn = shSyn(shrinkShVec(shcIn));
      [fRot,shcOut]= movePole(fIn,theta, false);
      
      for n=0:p
        row = n^2+1:(n+1)^2;
        R{n+1} = shcOut(row,p+1-n:p+1+n);
      end
      
      rReal =[]; rImag = [];
      for n=0:p
        rReal = [rReal;real(R{n+1}(:))];
        rImag = [rImag;imag(R{n+1}(:))];
      end
      writeData(fileName,num2str(nt, pstr),[rReal;rImag]);
    end
    Rall{nt} = R;
  end
  %printMsg('* The fast rotation data for p=%d were generate/read from file.\n',p);
end

function testKernelD()

  msg = {'For the sphere, normal vector is the eigenvector of','\n  the double layer with eigenvalue 1'};
  fprintf(['  ' cell2mat(msg) '\n ' repmat('-',1,length(msg{1})+2) '\n'])

  for p = [8 12 16]
    S = SurfaceSph(shape_gallery(p,''));
    Df = kernelD(S.geoProp.nor, S);
    %- For the sphere normal vector is the eigenvector of the double
    %layer with eigenvalue 1
    evalue = 1;
    err1 = dot(Df,S.geoProp.nor);
    err2 = norm(Df - evalue*err1.*S.geoProp.nor);
    fprintf('\n   error: %2.4e,\t%2.4e\n\n',max(abs(err1-evalue)),max(abs(err2)));
  end

  shape = 'tri_ap';
  %-Translations are eigenfunctions of the double layer with eigenvalue -1
  msg = 'Translations are eigenfunctions of the double layer with eigenvalue -1';
  fprintf(['  ' msg '\n ' repmat('-',1,length(msg)+2) '\n'])
  for p = [8 12 16]
    evalue = -1;
    S = SurfaceSph(shape_gallery(p,shape));
    d = vec3d(repmat([1 0 0], length(S.cart.x),1));
    Dd = kernelD(d, S);
    err1 = dot(Dd,d);
    err2 = norm(d - evalue*err1.*d);
    fprintf('\n   x translation error: %2.4e,\t%2.4e\n',max(abs(err1-evalue)),max(abs(err2)));

    d = vec3d(repmat([0 1 0], length(S.cart.x),1));
    Dd = kernelD(d,S);
    err1 = dot(Dd,d);
    err2 = norm(d - evalue*err1.*d);
    fprintf('   y translation error: %2.4e,\t%2.4e\n',max(abs(err1-evalue)),max(abs(err2)));

    d = vec3d(repmat([0 0 1], length(S.cart.x),1));
    Dd = kernelD(d,S);
    err1 = dot(Dd,d);
    err2 = norm(d - evalue*err1.*d);
    fprintf('   z translation error: %2.4e,\t%2.4e\n\n',max(abs(err1-evalue)),max(abs(err2)));
  end

  %- Solid rotations are eigenfunctions of the double layer with eigenvalue -1
  msg = 'Solid rotations are eigenfunctions of the double layer with eigenvalue -1';
  fprintf(['  ' msg '\n ' repmat('-',1,length(msg)+2) '\n'])
  for p = [8 12 16]
    evalue = -1;
    S = SurfaceSph(shape_gallery(p,shape));
    d = vec3d(repmat([1 0 0], length(S.cart.x),1));
    d = cross(d,S.cart);
    d = d./norm(d);

    Dd = kernelD(d,S);
    err1 = dot(Dd,d);
    err2 = norm(Dd - err1.*d);
    fprintf('\n   x rotation error: %2.4e,\t%2.4e\n',max(abs(err1-evalue)),max(abs(err2)));

    d = vec3d(repmat([0 1 0], length(S.cart.x),1));
    d = cross(d,S.cart);
    d = d./norm(d);

    Dd = kernelD(d,S);
    err1 = dot(Dd,d);
    err2 = norm(Dd - err1.*d);
    fprintf('   y rotation error: %2.4e,\t%2.4e\n',max(abs(err1-evalue)),max(abs(err2)));

    d = vec3d(repmat([0 0 1], length(S.cart.x),1));
    d = cross(d,S.cart);
    d = d./norm(d);

    Dd = kernelD(d,S);
    err1 = dot(Dd,d);
    err2 = norm(Dd - err1.*d);
    fprintf('   z rotation error: %2.4e,\t%2.4e\n\n',max(abs(err1-evalue)),max(abs(err2)));
  end
end