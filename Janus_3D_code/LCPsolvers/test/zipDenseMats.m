function  zipDenseMats
[dirname, ~] = setPaths();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% old params
    % srcDir = fullfile(dirname, '../data/amphi.lattice.n_5.p_8.cDist_2.5');
% correctAx = true;
% MC = 600;
% ps = [3,4,5,6,7,8];
% tols = 10.0.^(-5:-1:-8);
%% new params 
root = '/projects/niru8088/Spheroidal3D-collisions/Janus_3D_code/test/../resultsForRecord/'
prefix = 'amphi.lcp.lattice.n_5.p_8.cDist_3.lcpSlvr_proxquasinewton.polyDisperseRatio_0.2'
srcDir = fullfile(root, [prefix '.denseMats']);
correctA = false;
Nt = 50;
ps = [2,3,4,5,6,7,8];
tols = 10.0.^(-3:-1:-6);
ps = sort(ps);
tols = sort(tols);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['zipping mats in ' srcDir])
load(fullfile(root, [prefix '.mat']), 'lcp_list','Fparams')
C = {lcp_list.C};
F = {lcp_list.F};
b = arrayfun(@(lcp) lcp.b, lcp_list, 'UniformOutput',false);
disp(['- Found Nt=' num2str(Nt) ', ps=' num2str(ps) ', tols=' num2str(tols)]);
A = cell(Nt, numel(ps), numel(tols));
dt = cell(Nt, numel(ps), numel(tols));
for i = 1:Nt
    for j = 1:numel(ps)
        p = ps(j);
        if correctA 
            T = getTCorrector(p, C{i}, Fparams,F{i});
        end
        for k = 1:numel(tols)
            tol = tols(k);
            fn = ['ix_' num2str(i) '.p_' num2str(p) '.tol_' num2str(tol) '.mat'];
            filePath = fullfile(srcDir, fn);
            if ~exist(filePath,'file')
                disp(['- ' fn ' DNE'])
                continue
            end
            disp(['- ' fn]);
            res = load(filePath);
            if isfield(res, 'dt') && sum(res.dt > 0) ~= numel(res.dt)
                disp('-- fn needs to still be completed')
                disp(['--- only ' num2str(sum(res.dt > 0)) '/' num2str(numel(res.dt)) ' rows finished.'])
                continue
            elseif isempty(res.A)
                disp('--- No Collision Problem for this time')
                continue;
            elseif ~isfield(res, 'dt')
                disp('--- No saved dt: only fieldnames')
                continue;
            end
            dt{i,j,k} = res.dt;
            if correctA
                Abad = res.A;
                Acor = T(Abad);
                A{i,j,k} = Acor;
                AbadSym = (Abad + Abad')/2;
                AcorSym = (Acor + Acor')/2;
                err = norm(Abad - Acor) / norm(Acor);
                symErrBad = norm(Abad - AbadSym) / norm(AbadSym);
                symErrCor = norm(Acor - AcorSym) / norm(AcorSym);
                fprintf('--- relErr Abad - Acor %.4g\n', err)
                fprintf('--- sym relErr Abad %.4g\n', symErrBad)
                fprintf('--- sym relErr Acor %.4g\n', symErrCor)
            else 
                A{i,j,k} = res.A;
            end
        end
    end
end
%%
disp(['Saving to ' [srcDir '.allMats.mat'] ])
save([srcDir '.allMats.mat'], 'Nt', 'ps', 'tols', 'A', 'b', 'dt');

end
function T = getTCorrector(p,Ct, Fparams,Fhat)
n3 = size(Ct,1);
np = 2*p.*(p+1); 
%% Build the velocities the way Corona did for testing
Sc = SurfaceSph(shape_gallery(p,''));
Xp = reshape(Sc.cart.to_array,[],3);
Xg = repmat(Xp,n3,1) + reshape(repmat(Ct',np,1),3,[])';
%% W has to do with the sum over the surface of the sphere
[~, gwt]=g_grid(p+1);
wt = pi/p*repmat(gwt', 2*p, 1)./sin(gl_grid(p));
wt = wt(:);
WT = Sc.geoProp.W; WT= WT.*wt;
Wg = repmat(WT(:),n3,1);
W=[];
for i=1:n3
    idx = (1:np)+np*(i-1); %#ok<BDSCI> 
    X = Xg(idx,:);
    w = Wg(idx); 
    sumW = sum(w); 
    tau1 = sum(w.*X(:,3).^2)+sum(w.*X(:,2).^2); 
    tau2 = sum(w.*X(:,1).^2)+sum(w.*X(:,3).^2); 
    tau3 = sum(w.*X(:,2).^2)+sum(w.*X(:,3).^2); 

    W = [W; (1/sumW); (1/sumW); (1/sumW); ...
        (1/tau1); (1/tau2); (1/tau3)]; %#ok<AGROW> 
end
W = diag(W);

[colevent,collist,mindst,mindstsh] = LOCAL_check_collision_sph(Ct, Fparams); %#ok<ASGLU> 
if isempty(collist)
   T = @(Abad) Abad;
   return
end
%% Build F
ip = collist(:,1); jp = collist(:,2);
numF = length(ip);
% Compute vectors and normal vectors for pairs
R = Ct(ip,:)-Ct(jp,:);       %Ci - Cj numF x 3 (NIC: vector between particle pair centers)
NR = sqrt(sum(R.*R,2));      %|Ci-Cj| numF x 1 (NIC: distance between particle pairs centers)
Rhat = repmat(1./NR,1,3).*R; %eij = (Ci - Cj)/|Ci-Cj| (NIC: unit vectors between particle pairs)

%% Build A
F = zeros(6*n3,numF); % NIC: F maps contact to direction of force applied to a particular particle
for k=1:numF
    indi = (1:3)+6*(ip(k)-1);
    indj = (1:3)+6*(jp(k)-1);
    F(indi,k) = Rhat(k,:);
    F(indj,k) = -Rhat(k,:);
end

for k=numF+1:numF
    indi = (1:3)+6*(ipsh(k-numF)-1);
    F(indi,k) = -Ct(ipsh(k-numF),:)./norm(Ct(ipsh(k-numF),:));
end
assert(norm(F -Fhat) / norm(Fhat) < 1e-8);
%% Define the final corrector
T = @(Abad) F'*W*(F' \ Abad);
end


%% set defaults
%if ~exist('srcDir', 'var') || isempty(srcDir)
%     srcDir = fullfile(dirname, '../data/amphi.special.n_3.p_8.cDist_2.3');
%end

%% Set all outputs
% disp(['- looking for all *.mat files in ' srcDir])
% out = ls(fullfile(srcDir, '*.mat'));
% out = split(out);
% MC = -1;
% ps = [];
% tols = [];
% for i = 1:numel(out)
%     fn = out{i};
%     if isempty(out)
%         continue
%     end
%     prts = split(fn, [srcDir '/']);
%     fn = prts{end};
%     disp(['-- ' fn])
%     prts = split(fn, '.');
%     for j = 1:numel(prts)
%         prt = prts{j};
%         if contains(prt,'ix')
%             prt = split(prt, "_");
%             ix = str2double(prt{end});
%             disp(['--- ix = ' prt{end}])
%             if ix > MC
%                 disp(['---- MC = ' prt{end}])
%                 MC = ix;
%             end
%         elseif contains(prt,'p')
%             prt = split(prt, "_");
%             p = str2double(prt{end});
%             disp(['--- p = ' prt{end}])
%             if ~any(ps == p)
%                 ps(end+1) = p; %#ok<*AGROW> 
%                 disp(['---- ps = ' num2str(ps)])
%             end
%          elseif contains(prt,'tol')
%             prt = split(prt, "_");
%             tol = str2double(prt{end});
%             disp(['--- tol = ' prt{end}])
%             if ~any(tols == tol)
%                 tols(end+1) = tol;
%                 disp(['---- tols = ' num2str(tols)])
%             end
%         end
%     end
% end