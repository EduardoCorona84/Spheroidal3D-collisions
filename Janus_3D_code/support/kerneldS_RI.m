function [DMat,DRI] = kerneldS_RI(S, DRI, directIn)   

  persistent lastSurf %matVecHandle lastSurf
  
  if(nargin==0), testKernelD(); return; end
   
  pDirect = 1000; %empirical 
  direct = max(S.p) < pDirect;
 
  if(nargin>3), direct = directIn; end
 
  %if(direct && ( any([S.dblLayerOpStale]) || any([S ~= lastSurf]) ) )
  if direct
    lastSurf = S;
    [DMat,DRI] = kerneldSMatrix(S,DRI);
    S.dblLayerOpStale = false;
    %matVecHandle = @(den) multivec(DMat,den);
  else
    %matVecHandle = @(den) kerneldSMatFree(den, S );
  end
  
  %Df = matVecHandle(f);
end
    
function v=multivec(A,b)
l = length(A);
v{1:l} = zeros(size(b{1})); 

for ii=1:l
    v{ii} = A{ii}*b{ii};
end
end

function [DMat,DRI] = kerneldSMatrix(S,DRI)
  persistent p Ywt A Wsph

  if ( nargin > 0 && ~isempty(S) )
    printMsg('* Generating the direct double-layer matrix for p=%d.\n',S.p);
   
    np = length(S.cart.x);
    if (isempty(p) || S.p~=p)
      [Ywt,A,Wsph] = getPersistents(S.p);
      p = S.p;
    end
  
    W0 = S.geoProp.W./Wsph;
    prm = zeros(1,3*np); 
    prm(1:3:3*np) = 1:np; prm(2:3:3*np) = np+1:2*np; prm(3:3:3*np)=2*np+1:3*np;
    DMat = cell(3,3); 
    [DMat{:}] = deal(zeros(np));
    X = reshape(S.cart.to_array,np,3); 
    
    if isempty(DRI)
        DRI = zeros(np); 
        update=false; 
        Nr = reshape(S.geoProp.nor.to_array,np,3); 
        nx = Nr(:,1)'; ny = Nr(:,2)'; nz = Nr(:,3)'; 
    else
        update=true; 
    end
    ind_ref = reshape((1:np),p+1,2*p);

    for k = 1:2*p
      ind = circshift(ind_ref,[0 k-1]);
      ind = ind(:);
      
      for j = 1:p+1
        %- Rotating the pole to (j,k)
        R = A{j}(ind,ind);
        RX  =R*X; 
        xx = RX(:,1)'; 
        yy = RX(:,2)'; 
        zz = RX(:,3)'; 
        
        if ~update
        W = ((R*W0).*Wsph)';
        end
        ind_pole = (k-1)*(p+1) + j;

        %- The distance vector
        rx = xx - X(ind_pole,1);
        ry = yy - X(ind_pole,2);
        rz = zz - X(ind_pole,3);
        inv_rho = 1./sqrt(rx.^2 + ry.^2 + rz.^2);

        %- Normalizing
        rx = rx.*inv_rho;
        ry = ry.*inv_rho;
        rz = rz.*inv_rho;

        if ~update
        %- Laplace kernel and the scalar part
        NdotR = (nx(ind_pole).*rx + ny(ind_pole).*ry + nz(ind_pole).*rz);   
        g = -Ywt.*W.*NdotR.*inv_rho.*inv_rho;  
        
        %- Double layer Matrix
        DRI(ind_pole,:) = g;
        else
            g = DRI(ind_pole,:); 
        end
        
        DMat{1,1}(ind_pole,:) = (g.*rx.*rx)*R;
        DMat{1,2}(ind_pole,:) = (g.*rx.*ry)*R;
        DMat{1,3}(ind_pole,:) = (g.*rx.*rz)*R;

        DMat{2,2}(ind_pole,:) = (g.*ry.*ry)*R;      
        DMat{2,3}(ind_pole,:) = (g.*ry.*rz)*R;

        DMat{3,3}(ind_pole,:) = (g.*rz.*rz)*R;
        %{
        nx = Ywt.*W.*nx(ind_pole).*inv_rho.*inv_rho; 
        ny = Ywt.*W.*ny(ind_pole).*inv_rho.*inv_rho; 
        nz = Ywt.*W.*nz(ind_pole).*inv_rho.*inv_rho;
        DMat{1,1}(ind_pole,:) = (g.*(1+rx.*rx)   -2*nx.*rx)*R;
        DMat{1,2}(ind_pole,:) = (g.*(rx.*ry)-nx.*ry-ny.*rx)*R;
        DMat{1,3}(ind_pole,:) = (g.*(rx.*rz)-nx.*rz-nz.*rx)*R;
        DMat{2,2}(ind_pole,:) = (g.*(1+ry.*ry)   -2*ny.*ry)*R;      
        DMat{2,3}(ind_pole,:) = (g.*(ry.*rz)-ny.*rz-nz.*ry)*R;
        DMat{3,3}(ind_pole,:) = (g.*(1+rz.*rz)-2*nz.*rz)*R;
        %}
      end
    end
    DMat{2,1} = DMat{1,2};
    DMat{3,1} = DMat{1,3};
    DMat{3,2} = DMat{2,3};

    DMat = cell2mat(DMat);
    
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
