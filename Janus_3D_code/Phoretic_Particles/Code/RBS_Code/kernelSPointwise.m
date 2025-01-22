function  Sfs = kernelSPointwise(sig, S)
% SIGMA is real domain smooth scalar functions over the sphere. S is of class
% vesicle
  
  persistent p np qwNorth qwSouth
  
  if(nargin==0), testThis(); return; end
  
  %- Setting/correcting the persistent data
  if(isempty(p) || p~=S.p)
    p = S.p;
    printMsg('* Generating the point-wise single-layer data for p=%d.\n',p);
    np = length(S.cart.x);
    qwNorth = (1/8/pi)*SingularWeights(p); qwNorth = qwNorth(:)'; 
    qwSouth = flipud(reshape(qwNorth,p+1,[])); qwSouth = qwSouth(:)';
  end
  
  X0 = reshape(S.cart.to_array,[],3);
  rotData = [X0 sig];
  phiIdx = 1:2*p;
  T = SurfaceSph(S.cart);
  r = vec3d;
  
  for lat=1:floor(p/2)+1
    %Rotating to the desired lattitude
    rot = movePoleShcMatrix(rotData, lat, phiIdx);  
    
    for long=1:2*p
      T.cart = vec3d(rot{long}(:,1:3));
      Fs = T.geoProp.tensionOp(rot{long}(:,end));
      
      ll = (long-1)*(p+1);
      northIdx = ll + lat;
      southIdx = mod(np/2 + 1 + ll + p - lat, np) + 1;
      northPole = X0(northIdx,:);
      southPole = X0(southIdx,:);

      %- treating the north
      r.x = T.cart.x - northPole(1);
      r.y = T.cart.y - northPole(2);
      r.z = T.cart.z - northPole(3);

      invRho = norm(r);
      invRho = 1./invRho;
      
      den = Fs + invRho .* invRho .* dot(Fs,r).* r;
      den = T.geoProp.W .* invRho .* den;
      den = reshape(den.to_array,[], 3);
      Sfs(northIdx,:) = qwNorth * den;
      
      %- treating the south
      r.x = T.cart.x - southPole(1);
      r.y = T.cart.y - southPole(2);
      r.z = T.cart.z - southPole(3);

      invRho = norm(r);
      invRho = 1./invRho;
      
      den = Fs + invRho .* invRho .* dot(Fs,r).* r;
      den = T.geoProp.W .* invRho .* den;
      den = reshape(den.to_array,[], 3);
      Sfs(southIdx,:) = qwSouth * den;
    end
  end
  Sfs = vec3d(Sfs);
  
function testThis()
  
  p = 16;
  shape = 'ellipseZ';
  S = SurfaceSph(shape_gallery(p, shape));
  sig = S.geoProp.H;
  Fs = S.geoProp.tensionOp(sig);
  Sfs = kernelSPointwise(sig, S);
  SfRef = kernelSMatFree(Fs, S);
  SfRef = vec3d(SfRef);
  
  disp(max(norm(Sfs-SfRef)));

