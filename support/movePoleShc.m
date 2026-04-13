function sRot = movePoleShc(shc,thetaTarg,lambdaTarg) 
  
  if(nargin == 0), testMovePoleShc(); return;end
  
  shc = expandShVec(shc);
  [d1, d2]= size(shc);
  N = (sqrt(8*d1+1)-3)/4;
    
  %% lambdaTarg rotation
  mMat = exp(1i*repmat(-N:N,N+1,1)*lambdaTarg);
  shc = repmat(mMat(:),1,d2).*shc;

  %% thetaTarg rotation
  %Parametrization of the equator in the rotated from. We need 2*N+2
  %equispaced points in [0,2pi).
  phiR = (0:2*N+1)*2*pi/(2*N+2) + pi/(2*N+2);
  %Finding the corresponding material point
  [phi, theta] = cart2sph(cos(thetaTarg)*cos(phiR),sin(phiR),-sin(thetaTarg)*cos(phiR));
  %Mapping phi to [0,2pi) 
  phi = mod(phi,2*pi)';
  theta = pi/2-theta'; 
  %Coefficients for G
  c1 = sin(thetaTarg)*cos(theta).*cos(phi) - cos(thetaTarg)*sin(theta);
  c2 = sin(thetaTarg)*sin(phi)./sin(theta);
  
  theta = [theta(:);pi/2];
  phi = [phi(:);0];
  sRot = zeros(N+1,d2*(2*N+1));
  
  %Calculating f and g and finding their Fourier coefficients.
  for n=0:N    
    [Yn, Hn] = Ynm(n,[],theta,phi);
    
    P = ones(1, 2*N+1); P(N+1+(-n:n)) = Yn(end,:).'; Yn = Yn(1:end-1,:);
    Q = ones(1, 2*N+1); Q(N+1+(-n:n)) = Hn(end,:).'; Hn = Hn(1:end-1,:);

    Hn = -repmat(c1,1,2*n+1).*Hn + (c2*(1i*(-n:n))).*Yn;
    
    for idx=1:d2
      coeff = repmat(shc((N+(-n:n))*(N+1)+n+1,idx).',2*N+2,1);
      f = sum(coeff.*Yn,2);
      g = sum(coeff.*Hn,2);

      ff = fftshift(fft(f,[],1))/(2*N+2); ff = ff(2:end)';
      gg = fftshift(fft(g,[],1))/(2*N+2); gg = gg(2:end)';
      
      sRot(n+1, (idx-1) * (2*N+1) + (1:2*N+1) ) = (ff.*P + gg.*Q)./(P.^2+Q.^2);
    end
  end
  sRot = reshape(sRot,[],d2);
  sRot = shrinkShVec(sRot);

  function testMovePoleShc()
  np = 16;
  
  [u, v] = gl_grid(np); u = pi/2-u;
  rho = @(theta,phi) 1+.1*cos(6*theta);%.5*real(Ynm(4,3,theta,phi));
  [x, y, z] = sph2cart(v,u,rho(u,v));
  X =[x y z];
  theta = 3*pi/4; phi = pi/4;

  shc = shAna(X);
  
  XX = shSyn(movePoleShc(shc,theta,phi), true);
    
  subplot(1,2,1); plotb(X(:));
  subplot(1,2,2); plotb(XX(:));
  
  [phi, theta, r]=cart2sph(XX(:,1),XX(:,2),XX(:,3));
  e = abs(r - rho(theta,phi));
  disp(['error = ' num2str(max(e(:)))]);
  
 