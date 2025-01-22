function  Df = LkernelDMatFreeVec(I,J,S)
% KERNELDMATFREE - calculates the double-layer Laplace kernel matvec of the
% surface X, with area element size W, normal vector N, and density f.
% 
% EXAMPLE: 
%     p = 12;
%     S = boundary(p,'neck');
%     f = rand(3*2*p*(p+1),1);
%     Df = kernelSMatFree(f,S.cart.vecForm,S.geoProp.nor.vecForm, S.geoProp.W);
%
% The transformation matrices are built and stored as persistent for later
% use, therefore, the first call may be much slower than the later calls.
%

  persistent p np uRotMat vRotMat Wsph qw qwIm

  if(nargin==0), testDMF(); return; end

  dot3 = @(X,Y) X(1:end/3,:).*Y(1:end/3,:) + ...
         X(end/3+1:2*end/3,:).*Y(end/3+1:2*end/3,:) + ...
         X(2*end/3+1:end,:).*Y(2*end/3+1:end,:);

  % Number of columns to compute
  nf = size(I,1); 
  nCol = 8*nf;
  % Local (1 to np) and vesicle (1 to nv) indices
  IP = I(:,1); IV = I(:,2); I = I(:,3); 
  JP = J(:,1); JV = J(:,2); J = J(:,3); 
  
  % Build persistent rotation matrices and weights if they don't already
  % exist. 
  if(isempty(p) || p~=S(1).p)
    p = S(1).p;
    np = 2*p*(p+1);
    
    printMsg('* Generating the fast double-layer data for p=%d.\n',p);
    [u, v] = gl_grid(p);
    %Sperical domain area element
    Wsph = sin(u);
    v = reshape(v,p+1,[]);
    v = repmat(v(1,:),nCol,1); v = v(:)';  
    %Rotation matrix in the longitude direction
    vRotMat = 1i*repmat(-p:p,p+1,1);
    vRotMat = exp(shrinkShVec(vRotMat(:))*v);
    %Rotation matrix in the lattitude direction
    uRotMat = getMat(p);
    %Singular kernel integration weights
    % Here I changed the coefficient to 1/4/pi (Laplace)
    qw = (1/4/pi)*SingularWeights(p); qw = qw(:); 
    qwIm = flipud(reshape(qw,p+1,[])); qwIm = qwIm(:);
  else
    [~, v] = gl_grid(p); 
    v = reshape(v,p+1,[]);
    v = repmat(v(1,:),nCol,1); v = v(:)';  
    vRotMat = 1i*repmat(-p:p,p+1,1);
    vRotMat = exp(shrinkShVec(vRotMat(:))*v);
  end

  %- treating surfaces one by one, can be improved 
  X = reshape(S.cart.to_array,[],3);
  N = reshape(S.geoProp.nor.to_array,[],3);
  W = S.geoProp.W;
  
  XI = X(IT,:);    
  XJ = zeros(np,3*nf); NJ = XJ; WJ = zeros(np,nf); fJ = WJ; 
  
  for i=1:nf
     XJ(:,(3*(i-1)+1):3*i) = X((1:np)+np*(JV-1),:);
     NJ(:,(3*(i-1)+1):3*i) = N((1:np)+np*(JV-1),:);
     WJ(:,i) = W((1:np)+np*(JV-1));
     fJ(JP(i),i) = 1;  
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Here we have to build XI and XJ, WJ, NJ and fJ

  %Taking the spherical harmonic transform (f is now scalar)
  shc = shAna([XJ WJ./repmat(Wsph,1,nf) NJ fJ]);
  %Making copies for each v=const and aligning the pole with all of these.
  shc= vRotMat.*repmat(shc,1,2*p);
  %X,W,f index in the rotated matrix
  ii = reshape(1:nCol*2*p,nCol,[]);
  xInd = ii(1:3*nf,:); xInd = xInd(:);
  wInd = ii(3*nf+1:4*nf,:); wInd = wInd(:);
  fInd = ii(7*nf+1:8*nf,:); fInd = fInd(:);   
  nInd = ii(4*nf+1:7*nf,:); nInd = nInd(:);
    
  ii = reshape(1:np,p+1,[]);
  northHem = ii(1:floor(p/2)+1  ,:); 
  southHem = ii(floor(p/2)+2:end,:); 
  southHem = flipud(circshift(southHem,[0 -p]));

    Df = zeros(size(f));   
    for ii=1:floor(p/2)+1
      %Rotating to the desired lattitude
      shcRot = zeros((p+1)^2,size(shc,2)); 
      for jj = 0:p
        ind = jj^2+1:(jj+1)^2;
        shcRot(ind,:) = uRotMat{ii}{jj+1}*shc(ind,:);
      end
      %Finding the real-space values
      rot = shSyn(shcRot, false);
      %XP is the pole positions (north pole)
      XP = X(northHem(ii,:),:);
      XP = XP'; XP = XP(:)';
      %The distance in the physical space
      r = rot(:,xInd);
      r = r - repmat(XP,np,1);
      r = reshape(r,3*np,2*p);
      %Rotated surface element times the integration weights
      Wr = repmat(rot(:,wInd).*repmat(qw.*Wsph,1,2*p),nf,1);  
      %Rotated density
      fr = rot(:,fInd);
      fr = reshape(fr,nf*np,2*p);   
      %Rotated normal
      nr = rot(:,nInd);
      nr = reshape(nr,3*np,2*p);  
      %The distance and stokes kernel 
      rho = repmat(sqrt(dot3(r,r)),nf,1);
      rdotn = repmat(dot3(r,nr),nf,1); 
      rho = Wr.*fr.*rdotn./rho.^3;      
      %dumping the results   
      Df(northHem(ii,:),:) = reshape(sum(reshape(rho,np,[])).',nf,[]).';         
      
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
        Wr = repmat(rot(:,wInd).*repmat(qwIm.*Wsph,1,2*p),nf,1);  
        %The distance and stokes kernel
        rho = repmat(sqrt(dot3(r,r)),nf,1);
        rdotn = repmat(dot3(r,nr),nf,1); 
        rho = Wr.*fr.*rdotn./rho.^3;
        %dumping the results        
        Df(southHem(ii,:),:) = reshape(sum(reshape(rho,np,[])).',nf,[]).';           
      end
    
    %Df = Df(:);   
    Df = reshape(Df,[],nf);  
    
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
      printMsg('* Generating the fast rotation data for p=%d and target=%d.\n',p,nt);

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
  printMsg('* The fast rotation data for p=%d were generate/read from file.\n',p);

function testDMF()

  clear functions;
  p = 16;
  S = SurfaceSph(shape_gallery(p,'butternut'));
  f = vec3d( rand(3*2*p*(p+1),1));
  direct = true;
  opts = struct();
  F1 = S.doubleLayerMatVec(vec3d(f), opts, direct);
  F2 = kernelDMatFree(f,S);
  disp(norm(norm(F1 - F2) ) );
