function Sf = kernelS(f, S, opts, directIn)
% performs the singular single-layer Stokes. It call the correct function
% according to the size of S, spherical harmonics frequency, and options.

persistent GmatVec lastSurf

if(nargin==0), timeKernelS(); return; end

pDirect = 10000; %empirical form timeKernelS()
% Form the matrix only when there is one vesicle and p < pdirect
direct = max(S.p) < pDirect;

if( nargin< 3 ), opts = struct();end
if( nargin> 3 ), direct = directIn; end

if( isfield(opts,'mexUse') && opts.mexUse )
    error( 'The underlying mex code for singular Stokes is obsolete!' );
end

if direct
%if( direct && ( any([S.stokesOpStale]) || any(S ~= lastSurf) ) ) %%length(S) == 1 via direct
    lastSurf = S;
    for ii=1:length(S)
        GS{ii} = stokesKernelMat(to_array(S(ii).cart), S(ii).geoProp.W, S(ii).p);
        S(ii).stokesOpStale = false;
    end
    GmatVec = @(f) multivec(GS,f);
else
    GmatVec = @(f) kernelSMatFree(f, S);
end

%if(~isempty(f))
    Sf = GmatVec(f);
%else
%    Sf = [];
%end
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

function GS = stokesKernelMat(X, W, p)

persistent Ywt A Wsph np

printMsg('* Generating the direct Stokes matrix for p=%d.\n',p);
X = reshape(X,[],3);

if (isempty(Ywt) | size(X,1)~=np)
    printMsg('* Generating the direct Stokes integration data for p=%d.\n',p);
    Ywt = (1/8/pi)*SingularWeights(p);
    u = gl_grid(p);
    Wsph = sin(u);
    A = getDirectRotMats(p);
    np = 2*p*(p+1);
end

x = X(:,1); y  = X(:,2); z = X(:,3);
W0 = W./Wsph;
GS = cell(3,3); [GS{:}] = deal(zeros(numel(x), numel(x)));
gS = cell(3,1);

ind0 = reshape((1:numel(x))',p+1,2*p);
for k = 1:2*p
    ind = circshift(ind0,[0 k-1]); ind = ind(:);
    for j = 1:p+1
        R = A{j}(ind,ind);
        
        xx = (R*x)'; yy = (R*y)'; zz = (R*z)';
        
        W = ((R*W0).*Wsph)';
        indg = j+(p+1)*(k-1);
        
        inv_rho = 1./sqrt((xx - x(indg)).^2 + (yy - y(indg)).^2 + (zz - z(indg)).^2);
        %- Laplace kernel
        g = (Ywt.*W).*inv_rho;
        
        %- rest of Stokes kernel
        gS{1} = ((xx - x(indg)).*inv_rho);
        gS{2} = ((yy - y(indg)).*inv_rho);
        gS{3} = ((zz - z(indg)).*inv_rho);
        
        %- Stokes Matrix
        GS{1,1}(indg,:) = ( g.*(1 + gS{1}.^2) )*R;
        GS{2,2}(indg,:) = ( g.*(1 + gS{2}.^2) )*R;
        GS{3,3}(indg,:) = ( g.*(1 + gS{3}.^2) )*R;
        
        GS{1,2}(indg,:) = (g.*gS{1}.*gS{2})*R;
        GS{1,3}(indg,:) = (g.*gS{1}.*gS{3})*R;
        GS{2,3}(indg,:) = (g.*gS{2}.*gS{3})*R;
    end
end

GS{2,1} = GS{1,2}; GS{3,1} = GS{1,3}; GS{3,2} = GS{2,3};
GS = cell2mat(GS);

prm = zeros(1,3*np); 
prm(1:3:3*np) = 1:np; prm(2:3:3*np) = np+1:2*np; prm(3:3:3*np)=2*np+1:3*np;
GS = GS(prm,prm); 
end

function R = getDirectRotMats(p)

R = cell(p+1,1);
np = 2*p*(p+1);

for idx = 1:p+1
    fname = ['RotMat' num2str(idx) '-'];
    R{idx} = readData(fname, p, [np np]);
    
    if(isempty(R{idx}))        
        printMsg('* Generating direct rotation matrices for p=%d,idx=%d.\n',p,idx-1);
        u = gl_grid(p); R{idx} = movePole(eye(np),u(idx),0);
        writeData(fname, p, R{idx});
    end
end
printMsg('* Direct rotation matrices for p=%d were generated/read from file.\n',p);
end

function timeKernelS()
P = [6 12 16 24 32];
direct = [false true];
jmax = 15; %assumed number of usage
opts = struct('mexUse',false);

for kk=1:length(direct)
    for ii=1:length(P)
        printMsg([' ' repmat('-',1,75) '\n']);
        S = SurfaceSph(shape_gallery(P(ii),'square_ap'));
        
        %%To make sure that the data-file exist, and if does not exist it
        %is generated
        u = S.stokesMatVec(S.cart, opts, direct(kk));
        S.cart = S.cart;
        
        tic;
        u = S.stokesMatVec(S.cart, opts, direct(kk));
        setup{kk}(ii) = toc;
        disp(['  Setup : ' num2str(setup{kk}(ii))]);
        tic
        for jj=1:jmax
            u = S.stokesMatVec(S.cart, opts, direct(kk));
        end
        eval{kk}(ii) = toc;
        disp(['  Evaluation : ' num2str(eval{kk}(ii))]);
    end
end

subplot(1,3,1);
plot(P,setup{1},'o--k',P,setup{2},'.-k');
legend('Setup fast','Setup direct');

subplot(1,3,2);
plot(P,eval{1},'o--k',P,eval{2},'.-k');
legend('Eval fast','Eval direct');

subplot(1,3,3);
plot(P,eval{1} + setup{1},'o--k',P,eval{2} + setup{2},'.-k');
legend('Fast','Direct');
end


