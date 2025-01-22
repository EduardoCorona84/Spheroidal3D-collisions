function pot = DoubleAlltoAll(src, nor, den,trg)
% SRC, DEN, and NOR are ordered per vesicle, TRG is a single column vector
% assumed to be in ythe [x;y;z] format.

  if(nargin==0), testThis(); return; end
  dot3 = @(x,y) sum(reshape(x,[],3) .* reshape(y,[],3),2);
  
  srcIsTrg = false;
  if(nargin<4 || isempty(trg)) 
    nv = size(src, 2);
    trg = src'; 
    srcIsTrg = true;
  end
  trg = reshape(trg,[],3);
  src = src'; src = reshape(src,[],3);
  den = den'; den = reshape(den,[],3);
  nor = nor'; nor = reshape(nor,[],3);
  
  nsrc = size(src,1);
  ntrg = size(trg,1);
  pot = zeros(ntrg,3);
  r = zeros(nsrc,3);
  
  for ii=1:ntrg
    r(:,1) = src(:,1) - trg(ii,1);
    r(:,2) = src(:,2) - trg(ii,2);
    r(:,3) = src(:,3) - trg(ii,3);
        
    invRho = 1./sqrt(dot3(r,r));
    if(srcIsTrg)
      invRho(ii) = 0;
    end
    invRho = invRho.^5;
    rfn = dot3(r, den).*dot3(r, nor);
    rfn = rfn .* invRho;
    r   = repmat(rfn,1,3).*r;
    pot(ii,:) = sum(r,1);
  end

  if(srcIsTrg)
    pot = reshape(pot, nv, []);
    pot = pot';
  else
    pot = pot(:);
  end
  
function testThis()
  p = 12;
  S = SurfaceSph(shape_gallery(p, ''));
  src = S.cart.to_array;
  [~,gwt]=g_grid(p+1);
  wt = pi/p*repmat(gwt', 2*p, 1)./sin(gl_grid(p));
  
  den = repmat(S.geoProp.W .* wt, 3,1).*S.geoProp.nor.to_array;
  nor = S.geoProp.nor.to_array;
  trg = S.centerOfMass().to_array;
  pot = DoubleAlltoAll(src, den, nor, trg);
  disp(['  Integral over sphere: ' num2str(reshape(pot,[],3))]);
  
  nsrc = 10000;
  src = zeros(nsrc, 3);
  nor = zeros(nsrc, 3);
  den = zeros(nsrc, 3);
  src(:,1) = linspace(0,1, nsrc);
  nor(:,1) = 1;
  den(:,1) = 1/nsrc;
  src = src(:);
  den = den(:);
  nor = nor(:);
  trg = [2 0 0]';
  pot = DoubleAlltoAll(src, den, nor, trg);
  pot(1) = pot(1) + 1/2;
  disp(['  Error of the integral over a line: ' num2str(max(abs(pot)))]);

