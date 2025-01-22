function shc = shSynT(f)
% SHANA(F) - Calculate the spherical harmonics transform of the function
% set F. Each column of F is a function defined on the parameter domain
% defined by parDomain.
%
% SEE ALSO: PARDOMAIN, GRULE, YNM
%

  persistent LT LT_Freqs

  mem = whos('LT');
  if(mem.bytes > 4e8), LT = []; LT_Freqs = []; end

  %-- Extracting the size
  if(nargin==0), testShAna(); return;end
  [d1,d2] = size(f);
  p = (sqrt(2*d1+1)-1)/2;

  if(p~=fix(p))
    error(['The input should be defined on a Gauss-Legendre--uniform'...
           ' grid. Check the size of input.']);
  end

  %-- Allocating memory
  shc = zeros((2*p+1)*(p+1),d2);

  %-- Reshaping to match the parametrization grid structure
  f = reshape(f,p+1,2*p,d2);

  %-- DFT
  f = fft(f,[],2);
  f = circshift(reshape(f,2*p*(p+1),d2),p*(p+1));
  f(1:p+1,:) = 2*f(1:p+1,:);
  f = [f; f(1:p+1,:)];        

  %-- Legendre transform
  idx = find(LT_Freqs == p);

  if(isempty(idx))
    LT{end+1} = getLegMat(p);
    LT_Freqs(end+1) = p;
    idx = length(LT);
  end

  for m = -p:p-1 
    ind = (m+p)*(p+1)+1;
    shc(ind:ind+p,:) = LT{idx}{p+1+m}*f(ind:ind+p,:);
  end
  shc = shrinkShVec(shc);

function LT = getLegMat(p)

  printMsg('* Reading the Legendre Transform matrix for p=%d: ',p);
  [LT,rFlag] = readData('legTransT',p,[p+1 (2*p+1)*(p+1)]);

  if(~rFlag)
    printMsg('Failed\n');
    printMsg('* Generating the Legendre Transform matrix for p=%d\n',p);

    [el,wt]=g_grid(p+1);
    u = acos(el(:));

    LT  = cell(1,2*p+1);
    for n = 0:p
      Yn = Ynm(n, [], u, 0*u); %at v=0, Ynm = Pnm;
      for m = -n:n
        LT{p+1+m}(n+1,:) = Yn(:,n+1+m).';
      end
    end
    writeData('legTransT',p, cell2mat(LT));
    printMsg('* Stored generated LT matrices for p=%d\n', p);
  else
    printMsg('Successful\n');
    LT = mat2cell(LT, p+1, repmat(p+1,1,2*p+1));
  end

function testShAna()
%-- Self-test function for SHANA
  p = 8;
  [el,az] = gl_grid(p);

  for n=0:p
    for m=-n:n
      c = Ynm(n,m,el,az);
      shc = reshape(expandShVec(shAna(c(:))),p+1,2*p+1);

      hold on; cla;
      spy(shc>1e-8,'r');
      plot(p+m+1,n+1,'o');
      axis([1 2*p+1 1 p+1]);
      fprintf('%2d\t %+2d\t%5.5e\n',n, m, shc(n+1,p+m+1));
      hold off;
      pause(.1);
    end
  end

  rho = @(theta,phi) 1 +.5*exp(10*real(Ynm(2,1,theta,phi)).^2);
  x = sph2cart(az,el,rho(el,az));

  x1 = expandShVec(shAna(x));
  x2 = shAnaDirect(x,az,el);

  disp('Comaprison with the direct method:');
  disp(max(max(abs(x1-x2))));

function SHC = shAnaDirect(F,az,el)

  [d1, d2] = size(F);
  p = (sqrt(2*d1+1)-1)/2;

  shc = zeros(p+1,2*p+1);
  SHC = zeros((p+1)*(2*p+1),d2);

  for idx = 1:d2
    f = reshape(F(:,idx),p+1,2*p);
    for n = 0:p
      Yn = Ynm(n, [], el, az);
      for m = -n:n
        shc(n+1,p+1+m) = smoothQuad(f.*reshape(Yn(:,n+m+1),p+1,[]));
      end
    end
    SHC(:,idx) = shc(:);
  end

function val = smoothQuad(f)

  persistent wt
  n = size(f,1)-1;
  if(isempty(wt))
    [~, wt]=g_grid(n+1);
    wt = pi/n*wt;
  end
  val = wt*sum(f,2);
