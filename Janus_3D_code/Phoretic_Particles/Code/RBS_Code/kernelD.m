function Df = kernelD(f, S, ~, directIn)   

  persistent matVecHandle lastSurf
  
  if(nargin==0), testKernelD(); return; end
   
  pDirect = 1000; %empirical 
  direct = max(S.p) < pDirect;
 
  if(nargin>3), direct = directIn; end
 
  %if(direct && ( any([S.dblLayerOpStale]) || any([S ~= lastSurf]) ) )
  if direct
    lastSurf = S;
    for ii=1:length(S)
        DMat{ii} = kernelDMatrix(S(ii));
        S(ii).dblLayerOpStale = false;
    end
    matVecHandle = @(den) multivec(DMat,den);
  else
    matVecHandle = @(den) kernelDMatFree(den, S );
  end
  
  Df = matVecHandle(f);
end
    
function v=multivec(A,b)

if isempty(b)
    v = A{1}; 
else
l = length(A);
v{1:l} = zeros(size(b{1}));

for ii=1:l
    v{ii} = A{ii}*b{ii};
end
end
end

function DMat = kernelDMatrix(S)
  persistent p Ywt A Wsph

  if ( nargin > 0 && ~isempty(S) )
    printMsg('* Generating the direct double-layer matrix for p=%d.\n',S.p);
   
    np = length(S.cart.x);
    if (isempty(p) || S.p~=p)
      [Ywt,A,Wsph] = getPersistents(S.p);
      p = S.p;
    end
  
    W0 = S.geoProp.W./Wsph;
    DMat = cell(3,3); 
    [DMat{:}] = deal(zeros(np));
    ind_ref = reshape((1:np),p+1,2*p);

    for k = 1:2*p
      ind = circshift(ind_ref,[0 k-1]);
      ind = ind(:);
      
      for j = 1:p+1
        %- Rotating the pole to (j,k)
        R = A{j}(ind,ind);
        xx = (R*S.cart.x)'; nx = (R*S.geoProp.nor.x)';
        yy = (R*S.cart.y)'; ny = (R*S.geoProp.nor.y)';
        zz = (R*S.cart.z)'; nz = (R*S.geoProp.nor.z)';
        
        W = ((R*W0).*Wsph)';
        ind_pole = (k-1)*(p+1) + j;

        %- The distance vector
        rx = xx - S.cart.x(ind_pole);
        ry = yy - S.cart.y(ind_pole);
        rz = zz - S.cart.z(ind_pole);
        inv_rho = 1./sqrt(rx.^2 + ry.^2 + rz.^2);

        %- Normalizing
        rx = rx.*inv_rho;
        ry = ry.*inv_rho;
        rz = rz.*inv_rho;

        %- Laplace kernel and the scalar part
        g = (nx.*rx + ny.*ry + nz.*rz);
        g = g.*Ywt.*W.*inv_rho.*inv_rho;

        
        %- Double layer Matrix
        DMat{1,1}(ind_pole,:) = (g.*rx.*rx)*R;
        DMat{1,2}(ind_pole,:) = (g.*rx.*ry)*R;
        DMat{1,3}(ind_pole,:) = (g.*rx.*rz)*R;

        DMat{2,2}(ind_pole,:) = (g.*ry.*ry)*R;      
        DMat{2,3}(ind_pole,:) = (g.*ry.*rz)*R;

        DMat{3,3}(ind_pole,:) = (g.*rz.*rz)*R;
      end
    end
    DMat{2,1} = DMat{1,2};
    DMat{3,1} = DMat{1,3};
    DMat{3,2} = DMat{2,3};

    DMat = cell2mat(DMat);
    
    prm = zeros(1,3*np); 
    prm(1:3:3*np) = 1:np; prm(2:3:3*np) = np+1:2*np; prm(3:3:3*np)=2*np+1:3*np;
    DMat = DMat(prm,prm);
  end
end

function [Ywt,A,Wsph] = getPersistents(p)

  printMsg('* Generating the direct double-layer integration data for p=%d.\n',p);
  Ywt = -(3/4/pi)*SingularWeights(p);
  u = gl_grid(p);
  Wsph = sin(u);
  np = 2*p*(p+1);

  A = cell(p+1,1);
  for idx = 1:p+1
    fname = ['RotMat' num2str(idx) '-'];
    A{idx} = readData(fname, p, [np np]);
    
    if(isempty(A{idx}))
      printMsg('* Generating the direct rotation matrices for p=%d, idx=%d.\n',p,idx-1);
      A{idx} = movePole(eye(np),u(idx),0);
      writeData(fname, p, A{idx});
    end
  end
  printMsg('* Direct rotation matrices for p=%d were generated/read from file.\n',p);
end

function testKernelD()

  msg = {'For the sphere, normal vector is the eigenvector of','\n  the double layer with eigenvalue 1'};
  fprintf(['  ' cell2mat(msg) '\n ' repmat('-',1,length(msg{1})+2) '\n'])

  for p = [8 12 16]
    S = SurfaceSph(shape_gallery(p,''));
    Df = kernelD(S.geoProp.nor, S);
    %- For the sphere normal vector is the eigenvector of the double
    %layer with eigenvalue 1
    evalue = 1;
    err1 = dot(Df,S.geoProp.nor);
    err2 = norm(Df - evalue*err1.*S.geoProp.nor);
    fprintf('\n   error: %2.4e,\t%2.4e\n\n',max(abs(err1-evalue)),max(abs(err2)));
  end

  shape = 'tri_ap';
  %-Translations are eigenfunctions of the double layer with eigenvalue -1
  msg = 'Translations are eigenfunctions of the double layer with eigenvalue -1';
  fprintf(['  ' msg '\n ' repmat('-',1,length(msg)+2) '\n'])
  for p = [8 12 16]
    evalue = -1;
    S = SurfaceSph(shape_gallery(p,shape));
    d = vec3d(repmat([1 0 0], length(S.cart.x),1));
    Dd = kernelD(d, S);
    err1 = dot(Dd,d);
    err2 = norm(d - evalue*err1.*d);
    fprintf('\n   x translation error: %2.4e,\t%2.4e\n',max(abs(err1-evalue)),max(abs(err2)));

    d = vec3d(repmat([0 1 0], length(S.cart.x),1));
    Dd = kernelD(d,S);
    err1 = dot(Dd,d);
    err2 = norm(d - evalue*err1.*d);
    fprintf('   y translation error: %2.4e,\t%2.4e\n',max(abs(err1-evalue)),max(abs(err2)));

    d = vec3d(repmat([0 0 1], length(S.cart.x),1));
    Dd = kernelD(d,S);
    err1 = dot(Dd,d);
    err2 = norm(d - evalue*err1.*d);
    fprintf('   z translation error: %2.4e,\t%2.4e\n\n',max(abs(err1-evalue)),max(abs(err2)));
  end

  %- Solid rotations are eigenfunctions of the double layer with eigenvalue -1
  msg = 'Solid rotations are eigenfunctions of the double layer with eigenvalue -1';
  fprintf(['  ' msg '\n ' repmat('-',1,length(msg)+2) '\n'])
  for p = [8 12 16]
    evalue = -1;
    S = SurfaceSph(shape_gallery(p,shape));
    d = vec3d(repmat([1 0 0], length(S.cart.x),1));
    d = cross(d,S.cart);
    d = d./norm(d);

    Dd = kernelD(d,S);
    err1 = dot(Dd,d);
    err2 = norm(Dd - err1.*d);
    fprintf('\n   x rotation error: %2.4e,\t%2.4e\n',max(abs(err1-evalue)),max(abs(err2)));

    d = vec3d(repmat([0 1 0], length(S.cart.x),1));
    d = cross(d,S.cart);
    d = d./norm(d);

    Dd = kernelD(d,S);
    err1 = dot(Dd,d);
    err2 = norm(Dd - err1.*d);
    fprintf('   y rotation error: %2.4e,\t%2.4e\n',max(abs(err1-evalue)),max(abs(err2)));

    d = vec3d(repmat([0 0 1], length(S.cart.x),1));
    d = cross(d,S.cart);
    d = d./norm(d);

    Dd = kernelD(d,S);
    err1 = dot(Dd,d);
    err2 = norm(Dd - err1.*d);
    fprintf('   z rotation error: %2.4e,\t%2.4e\n\n',max(abs(err1-evalue)),max(abs(err2)));
  end
end