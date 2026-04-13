function  fRot = movePoleShcMatrix(f, thetaIdx, phiIdx, isReal)
% f is a matrix such that each column is a function over the sphere on a
% Gauss-uniform grid returned by PARDOMAIN(). THETAIDX and PHIIDX are vectors
% containing the index of the desired target points. The pole of the input
% function is rotated to the tensor product of THETAIDX and PHIIDX, i.e. all
% pairs of combination of the members of THETAIDX and PHIIDX. At return, fRot
% holds a cell of size LENGTH(THETAIDX) X LENGTH(PHIIDX) and each element is the
% rotated matrix of f.
  
  persistent p uRotMat v vRotMat
  
  if(nargin<4), isReal = true; end
  np = size(f,1);
  pIn = (sqrt(2*np+1)-1)/2;
  
  %- Setting/correcting the persistent data
  if(isempty(p) || p~=pIn)
    p = pIn;
    printMsg('* Generating the SHC rotation data for p=%d.\n',p);
    [u,v] = gl_grid(p);
    v = reshape(v, p+1, []);
    v = v(1,:)';
    vRotMat = 1i*repmat(-p:p,p+1,1);
    vRotMat = shrinkShVec(vRotMat(:));
    uRotMat = getShRotMat(p);
  end

  nTheta = length(thetaIdx);
  nPhi = length(phiIdx);
  nFun = size(f,2);
  
  %- Generating the longitude rotations and allocating space
  vv = v(phiIdx); 
  vv = repmat(vv(:)', nFun, 1);
  vv = vv(:)';
  vMat = vRotMat * vv;
  vMat = exp(vMat);
  
  fRot = cell(nTheta, nPhi);
  
  %Taking the spherical harmonic transform
  shcRef = shAna(f);
    
  for ii = 1:nTheta
    fTmp = repmat(shcRef, 1, nPhi); 
    
    %Rotating to the desired longitude(s)
    fTmp = vMat .* fTmp;
    
    %Rotating to the latitude  
    for jj = 0:p
      ind = jj^2+1:(jj+1)^2;
      fTmp(ind,:) = uRotMat{thetaIdx(ii)}{jj+1}*fTmp(ind,:);
    end
    
    fRot(ii,:) = mat2cell(shSyn(fTmp, isReal), np, nFun*ones(1,nPhi));
  end
