
load('/projects/niru8088/Spheroidal3D-collisions/Janus_3D_code/data/amphi.lattice.n_5.p_8.cDist_2.5.mat')
load('/projects/niru8088/Spheroidal3D-collisions/Janus_3D_code/LCPsolvers/data/amphi.lattice.n_5.p_8.cDist_2.5.mat', 'lcp_list' ,'Fparams');
rd = Fparams.parbd.rd; 
n3 = length(rd); 
diam = Fparams.parbd.diam; %diam(i,j) = r_i + r_j
mxrd = Fparams.parbd.mxrd; %max(r_i,r_j)
eps = Fparams.parbd.eps;
Nb = Fparams.parbd.Nb; 
for i = 1:numel(lcp_list)
    if ~isempty(lcp_list(i).A)
        continue
    end
    C = Ct{i};
    [colevent,collist,mindst,mindstsh] = LOCAL_check_collision_sph(C,Fparams);
    if isempty(collist)
        continue
    end
    %% Setup (build A and b)
    ip = collist(:,1); jp = collist(:,2);
    numFS = 0;
    numF = length(ip);

    % Compute vectors and normal vectors for pairs
    R = C(ip,:)-C(jp,:);       %Ci - Cj numF x 3 (NIC: vector between particle pair centers)
    NR = sqrt(sum(R.*R,2));      %|Ci-Cj| numF x 1 (NIC: distance between particle pairs centers)
    Rhat = repmat(1./NR,1,3).*R; %eij = (Ci - Cj)/|Ci-Cj| (NIC: unit vectors between particle pairs)

    %% Build A
    F = zeros(6*n3,numF+numFS); % NIC: F maps contact to direction of force applied to a particular particle

    for k=1:numF
        indi = (1:3)+6*(ip(k)-1);
        indj = (1:3)+6*(jp(k)-1);
        F(indi,k) = Rhat(k,:);
        F(indj,k) = -Rhat(k,:);
    end

    for k=numF+1:numF+numFS
        indi = (1:3)+6*(ipsh(k-numF)-1);
        F(indi,k) = -C(ipsh(k-numF),:)./norm(C(ipsh(k-numF),:));
    end

    %% Build constant vector b:
    % Compute (1/dt)*phi
    phib = zeros(numF+numFS,1);
    if numF>0
        phib(1:numF) = (1/dt)*(NR-diam(ip+n3*(jp-1))-eps*mxrd(ip+n3*(jp-1)));
    end

    bvec = phib + real((F.')*VW{i}(:));
    lcp_list(i).b = bvec;
    lcp_list(i).F = F;
    lcp_list(i).C = C;
end
save('/projects/niru8088/Spheroidal3D-collisions/Janus_3D_code/LCPsolvers/data/amphi.lattice.n_5.p_8.cDist_2.5.mat', 'lcp_list' ,'Fparams');