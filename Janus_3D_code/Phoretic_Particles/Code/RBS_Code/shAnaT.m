function f = shAnaT(shc,isReal)
% F = SHSYN(SHC,isReal) - Calculates the inverse spherical transform for the
% given set of spherical harmonics coefficients SHC. Each column of SHC
% is the coefficients of a function. isReal is the flag that indicates
% that the desired function F is real valued.
%
% SEE ALSO: SHANA, YNM.
%

  persistent LTI LTI_Freqs

  mem = whos('LTI');
  if(mem.bytes > 4e8), LTI = []; LTI_Freqs = []; end

  %-- Checking the input
  if(nargin==0), testShSyn(); return; end
  if(nargin==1), isReal = false; end

  shc = expandShVec(shc);
  %-- Calculating the sizes
  [d1, d2]= size(shc);
  p = (sqrt(8*d1+1)-3)/4;

  if(p~=fix(p))
    error(['The input should be defined on a Gauss-Legendre--uniform'...
           ' grid. Check the size of input.']);
  end

  %-- Inverse Legendre transform
  idx = find(LTI_Freqs == p);

  if(isempty(idx))
    LTI{end+1} = getLegSynMat(p);
    LTI_Freqs(end+1) = p;
    idx = length(LTI);
  end

  f = zeros(2*p*(p+1),d2); ftmp = [f ; zeros(p+1,d2)]; 
  for m = -p:p %The last frequency is extra       
    ind = (m+p)*(p+1)+1;
    ftmp(ind:ind+p,:) = LTI{idx}{p+1+m}*shc(ind:ind+p,:);
  end

  f = ftmp(1:size(f,1),:);    
  f(1:p+1,:) = (ftmp(1:p+1,:)+ftmp(size(f,1)+1:end,:))/2;          
  f = circshift(f,-p*(p+1));
  f = reshape(f,p+1,2*p,d2);
  f = ifft(f,[],2)*2*p;
  if(isReal), f = real(f);end
  f = (pi/p)*reshape(f,[],d2);  

function LTI = getLegSynMat(p)

  printMsg('* Reading the inverse Legendre Transform matrix for p=%d: ',p);
  [LTI,  rFlag] = readData('invLegTransT',p,[p+1 (2*p+1)*(p+1)]);

  if(~rFlag)
    printMsg('Failed\n');
    printMsg('* Generating the inverse Legendre Transform matrix for p=%d\n',p);

    [el,wt] = g_grid(p+1);
    wt = wt.'; 
    u = acos(el(:));

    LTI  = cell(1,2*p+1);
    for n = 0:p
      Yn = Ynm(n, [], u, 0*u); %at v=0, Ynm = Pnm;
      for m = -n:n
        LTI{p+1+m}(:,n+1) = wt.*Yn(:,n+1+m); 
      end
    end
    writeData('invLegTransT',p, cell2mat(LTI));
    printMsg('* Stored generated ILT matrices for p=%d\n', p);
  else
    LTI = mat2cell(LTI, p+1, repmat(p+1,1,2*p+1));
    printMsg('Successful\n');
  end

function testShSyn()

  N = 8;
  [el, az] = gl_grid(N);

  sep = '-------------------------------\n';
  fprintf([ sep '   n        m     error\n' sep]);
  for n=0:N
    for m=-n:n
      f = Ynm(n,m,el,az);
      ff = shSyn(shAna(f));
      fprintf('   %d        %+d     %2.4e\n',n,m,max(max(abs(f-ff))));
    end
  end
  fprintf([sep '  A more general function \n' sep]);
  f = 1+cos(6*el);
  ff = shSyn(shAna(f),1);
  fprintf(['   error = %2.4e\n' sep],max(max(abs(f-ff))));
