function  Df = kernelDPointwise(f,X,N,W)
% KERNELDMATPOINTWISE - calculates the double-layer stokes kernel matvec of the
% surface X, with area element size W, normal vector N, and density f.
% 
% EXAMPLE: 
%     p = 12;
%     S = boundary(p,'neck');
%     f = rand(3*2*p*(p+1),1);
%     Df = kernelSMatPointwise(f,S.cart.vecForm,S.geoProp.nor.vecForm, S.geoProp.W);
%
% The transformation matrices are built and stored as persistent for later
% use, therefore, the first call may be much slower than the later calls.
%

  persistent p np ind uRotMat vRotMat v Wsph qw qwIm

  if(nargin==0), testDMF(); return; end
  dot3 = @(X,Y) X(1:end/3,:).*Y(1:end/3,:) + ...
         X(end/3+1:2*end/3,:).*Y(end/3+1:2*end/3,:) + ...
         X(2*end/3+1:end,:).*Y(2*end/3+1:end,:);

  X = reshape(X,[],3);
  f = reshape(f,[],3);
  N = reshape(N,[],3);
  
  nCol = 10;
  %Building persistent variables
  if(isempty(ind) || np~=size(X,1))
    np = size(X,1);
    p = (sqrt(2*np+1)-1)/2;
    printMsg('* Generating the point-wise double-layer data for p=%d.\n',p);
    %The matrices act on each moment of Ynm, that does not agree with our
    %ordering, ind is the reordering index.
    ind = reshape(1:(p+1)*(2*p+1),p+1,[]);
    ind = ind'; ind = ind(:);
    [u, v] = gl_grid(p);
    %Sperical domain area element
    Wsph = sin(u);
    v = reshape(v,p+1,[]);
    v = v(1,:);
    %Rotation matrix in the longitude direction
    vRotMat = 1i*repmat(-p:p,p+1,1);
    vRotMat = repmat(vRotMat(ind),1,nCol);
    %Rotation matrix in the lattitude direction
    uRotMat = getMat(p);
    %Singular kernel integration weights
    qw = -(3/4/pi)*SingularWeights(p); qw = qw(:); 
    qwIm = flipud(reshape(qw,p+1,[])); qwIm = qwIm(:);
  end

  %Taking the spherical harmonic transform
  shc0 = expandShVec(shAna([X W./Wsph f N]));
  shc0 = shc0(ind,:);
    
  Df = zeros(size(X));
  for longIdx=1:2*p
    %Aligning the meridian
    shc = exp(vRotMat*v(longIdx));
    shc = shc.*shc0;
    
    for ii=1:floor(p/2)+1
      %Rotating to the desired lattitude
      for jj = 1:p+1
        a = (jj-1)*(2*p+1)+1;
        shcRot(a:a+2*p,:) = uRotMat{ii}{jj}*shc(a:a+2*p,:);
      end
      %shuffling back to our standard order
      shcRot(ind,:) = shcRot; 
      
      %Finding the real-space values
      rot = shSyn(shcRot, false);

      %XN is the pole positions (north pole)
      ll = (longIdx-1)*(p+1);
      northPole = ll + ii;
      southPole = mod(np/2 + 1 + ll + p - ii, np) + 1;
      XN = X(northPole,:);
      XS = X(southPole,:);
      %The distance in the physical space
      r = rot(:,1:3);
      r = r - repmat(XN,np,1);
      r = r(:);     
      %Rotated surface element times the integration weights
      Wr = rot(:,4).*qw.*Wsph;
      Wr = repmat(Wr,3,1);
      %Rotated density
      fr = rot(:,5:7);
      fr = fr(:);
      %Rotated normal
      nr = rot(:,8:10);
      nr = nr(:);

      %The distance and stokes kernel
      rho = sqrt(dot3(r,r)); 
      rho = repmat(rho,3,1);
      rho = Wr.*(repmat(dot3(r,fr).*dot3(r,nr),3,1))./rho.^5;
      rho = rho.*r;
      %dumping the results
      Df(northPole,:) = sum(reshape(rho,[],3));
     
      %Treating the south pole at the same time
      if(ii~=p/2+1)
        %The distance in the physical space
        r = rot(:,1:3);
        r = r - repmat(XS,np,1);
        r = r(:);
        %Rotated surface element times the integration weights
      Wr = rot(:,4).*qwIm.*Wsph;
      Wr = repmat(Wr,3,1);
      %The distance and stokes kernel
      rho = sqrt(dot3(r,r));
      rho = repmat(rho,3,1);
      rho = Wr.*(repmat(dot3(r,fr).*dot3(r,nr),3,1))./rho.^5;
      rho = rho.*r;
      %dumping the results
      Df(southPole,:) = sum(reshape(rho,[],3));
      end
    end
  end
  Df = Df(:);

function Rall = getMat(p)

  nmat = floor(p/2)+1;
  Rall = cell(nmat,1);
  
  for nt=1:nmat 
    R = [];
    R = cell(p+1,1);
    fileName = ['rotMat_' num2str(p) '_'];
    matSize = (2*p+1)^2;
    count = (p+1)*matSize;

    RM = readData(fileName,nt,2*count);
    if(~isempty(RM))
      rReal = RM(1:count);
      rImag = RM(count+1:end);
      for n=0:p
        ind = n*matSize;
        R{n+1} = reshape(complex(rReal(ind+1:ind+matSize), ...
                                 rImag(ind+1:ind+matSize)),2*p+1,2*p+1);
      end
    else
      u = g_grid(p);
      theta = u(nt);
      printMsg('* Generating the fast rotation data for p=%d.\n',p);
      for n=0:p
        R{n+1} = zeros(2*p+1);
      end

      shcIn = zeros((p+1)*(2*p+1),2*p+1);
      for m=1:2*p+1
        ind = (m-1)*(p+1)+1;
        shcIn(ind:ind+p,m) = 1;
      end
      
      fIn = shSyn(shrinkShVec(shcIn), false);
      [~, shcOut]= movePole(fIn,theta, false);
      for m=-p:p
        temp = reshape(shcOut(:,p+1+m),p+1,[]);
        for n=0:p
          R{n+1}(:,p+m+1) = temp(n+1,:);
        end
      end

      rReal =[]; rImag = [];
      for n=0:p
        rReal = [rReal;real(R{n+1}(:))];
        rImag = [rImag;imag(R{n+1}(:))];
      end
      writeData(fileName,nt,[rReal;rImag]);
    end
    Rall{nt} = R;
  end
  printMsg('* The fast rotation data for p=%d were generate/read from file.\n',p);
  
function testDMF()

  %clear functions;
  p = 16;
  S = SurfaceSph(shape_gallery(p,'butternut'));
  f = rand(3*2*p*(p+1),1);
  direct = true;
  F1 = S.doubleLayerMatVec(vec3d(f),direct);
  F2 = kernelDPointwise(f, S.cart.to_array,S.geoProp.nor.to_array,S.geoProp.W);
  disp(max(norm(F1-vec3d(F2))));
