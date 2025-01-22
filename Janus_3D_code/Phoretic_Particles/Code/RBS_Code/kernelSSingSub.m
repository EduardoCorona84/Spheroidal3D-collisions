function Sf = kernelSSingSub(sigma, S)

  persistent wt p

  if (nargin==0), testKSSS(); return; end
  
  if (isempty(wt) || (~isempty(S) && p ~= S.p) )
    p = S.p;
    [trash gwt]=grule(p+1);
    wt = (1/8/pi) * pi/p*repmat(gwt', 2*p, 1)./sin(parDomain(p));
    wt = wt(:)';
  end
    
  f = S.geoProp.tensionOp(sigma); %%Non-normal vector does not work
  %f = sigma .*S.geoProp.H .* S.geoProp.nor;
  r = f; Sf = f;
  
  for ii=1:length(f.x)
    r.x = S.cart.x - S.cart.x(ii);
    r.y = S.cart.y - S.cart.y(ii);
    r.z = S.cart.z - S.cart.z(ii);
    
    invRho = norm(r);
    invRho = 1./invRho;
    invRho(ii) = 0;

    den = f + (invRho .* invRho .* dot(r, f)) .* r;
    den = S.geoProp.W .* invRho .* den;

    Sf.x(ii) = wt * den.x; 
    Sf.y(ii) = wt * den.y; 
    Sf.z(ii) = wt * den.z; 
  end
  
  
function testKSSS()

  shape = 'neck';
  for p = [12 16 24 32]
    clear functions
    S = boundary(p, shape);
  
    sigma = shSyn(rand((p+1)^2,1), 1); 
    Fs = sigma .* S.geoProp.H .*S.geoProp.nor;
  
    SfRef = S.stokesMatVec(Fs);
    Sf    = kernelSSingSub(sigma, S);
  
    disp(max(norm(Sf-SfRef))/max(norm(SfRef)));
  end
  
