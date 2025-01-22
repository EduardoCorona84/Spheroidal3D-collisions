function shcRot = movePoleFast(shc,beta) 
  
  if(nargin==0), testMovePoleFast(); return;end
  
  shc = expandShVec(shc);
  [d1, d2]= size(shc);
  p = (sqrt(8*d1+1)-3)/4;
   
  %- beta rotation Parametrization of the equator in the rotated from. We
  %- need 2*p+2 equi-spaced points in [0,2pi).
  phiR = (0:2*p+1)*2*pi/(2*p+2) + pi/(2*p+2);
  
  %- Finding the corresponding material point and mapping phi to [0,2pi) 
  [phi, theta] = cart2sph(cos(beta)*cos(phiR),sin(phiR),-sin(beta)*cos(phiR));
  theta = pi/2 - theta';
  phi = mod(phi,2*pi)';
  
  %- Coefficients for G
  c1 = sin(beta)*cos(theta).*cos(phi) - cos(beta)*sin(theta);
  c2 = sin(beta)*sin(phi)./sin(theta);
    
  %- Calculating F and G and finding their Fourier coefficients.
  shc = reshape(shc, p+1, d2*(2*p+1) );
  %- theta=pi/2 and phi = 0 gives the Legendre polynomials on the equator
  theta = [theta(:); pi/2];
  phi = [phi(:); 0];
  
  %- Allocating memory
  shcRot = zeros(d2*2*p*(2*p+1), p+1);
  
  for n=0:p    
    [Yn, Hn] = Ynm(n,[],theta(:),phi(:));
  
    %% Used to recover the rotated coefficients
    P0 = ones(1, 2*p+1); P0(p+1-n:p+1+n) = Yn(end,:).'; Yn = Yn(1:end-1,:);
    Q0 = ones(1, 2*p+1); Q0(p+1-n:p+1+n) = Hn(end,:).'; Hn = Hn(1:end-1,:);
    
    if ( n == p )
      Yn = Yn(:,2:end);
      Hn = Hn(:,2:end);
      Hn = -repmat(c1,1,2*n).*Hn + (c2*(1i*(-n+1:n))).*Yn;
    else
      Hn = -repmat(c1,1,2*n+1).*Hn + (c2*(1i*(-n:n))).*Yn;
    end
       
    for ii=0:d2-1
      %% Padding the coefficients to feed to 2p point fft
      F = zeros(2*p+2,2*p);
      G = zeros(2*p+2,2*p);
  
      if ( n == p )        
        coeff = repmat(shc(n+1, (2*p+1)*ii + (2:2*p+1) ),2*p+2,1);
                
        F = coeff.*Yn;
        G = coeff.*Hn;
      else
        coeff = repmat(shc(n+1,  (2*p+1)*ii+p+1-n:(2*p+1)*ii+p+1+n),2*p+2,1);
        
        F(:,p+1-n:p+1+n) = coeff.*Yn;
        G(:,p+1-n:p+1+n) = coeff.*Hn;
      end

      F = ifft(ifftshift(F, 2), [], 2) * 2 * p;
      G = ifft(ifftshift(G, 2), [], 2) * 2 * p;
      
      F = fftshift(fft(F,[],1),1)/(2*p+2);
      G = fftshift(fft(G,[],1),1)/(2*p+2);
      
      %%calculating the new coeffiecients
      F = F(2:end,:); F = F(:)';
      G = G(2:end,:); G = G(:)';
      P = repmat(P0,1, 2*p);
      Q = repmat(Q0,1, 2*p);
      
      shcRot(ii*2*p*(2*p+1)+1:(ii+1)*2*p*(2*p+1),n+1) = (F .* P + G .* Q) ./ (P.*P+Q.*Q);
    end
  end
  shcRot = reshape(shcRot',d1,[]);
  shcRot = shrinkShVec(shcRot);

  
function testMovePoleFast()
  
  p = 16;
  [u, v] = gl_grid(p); u = pi/2-u;
  rho = @(theta,phi) 1+.1*cos(6*theta);
  [x, y, z] = sph2cart(v,u,rho(u,v));
  X =[x y z];
  
  theta = pi/3;
  shc = shAna(X);
  
  shcRot = movePoleFast(shc, theta);
  XRot = reshape(shSyn(shcRot, true),[],3);
  
  np = 2*p*(p+1);
  for ii=1:2*p
    XX = XRot((ii-1)*np+1:ii*np,:);
    plotb(XX(:));
    drawnow;
    pause(.1);
    
    [phi, theta, r]=cart2sph(XX(:,1),XX(:,2),XX(:,3));
    e = abs(r - rho(theta,phi));
    disp(['error = ' num2str(max(e(:)))]);
  end
  