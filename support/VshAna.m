function Vshc = VshAna(F,type)
% VshAna(F) - Calculate the vector spherical harmonics transform of the 
% vector field set F. Each column of F is a vector field defined on the parameter domain
% defined by parDomain.
% Input is assumed to be of the form [Fx;Fy;Fz] where F is a vector field
% evaluated on a Gauss-Legendre uniform grid of order p
%
% SEE ALSO: PARDOMAIN, GRULE, YNM, VNM
%

Vshc=[]; 

% Extract size (6*p*(p+1) x d2)
if(nargin==0), testVShAna(); return;end
[d1,d2] = size(F);
np = d1/3;
p = (sqrt(2*np+1)-1)/2;
sp=(p+1)^2;

if(p~=fix(p))
    error(['The input should be defined on a Gauss-Legendre--uniform'...
           ' grid. Check the size of input.']);
end

% Setup 
[u,v]=gl_grid(p); 
% Normal Unit vector
er = [sin(u).*cos(v) sin(u).*sin(v) cos(u)];
% Tangential Unit vectors 
ev = [-sin(v) cos(v) zeros(np,1)]; 
eu = [cos(u).*cos(v) cos(u).*sin(v) -sin(u)];  

dotp = @(x,y) x(1:np,:).*conj(y(1:np,:)) + x(np+1:2*np,:).*conj(y(np+1:2*np,:)) + x(2*np+1:3*np,:).*conj(y(2*np+1:3*np,:));

%Separate F into normal, u and v components
Fr = dotp(F,repmat(er(:),1,d2)); 
Fv = dotp(F,repmat(ev(:),1,d2)); 
Fu = dotp(F,repmat(eu(:),1,d2)); 

%Analyze r component to obtain radial coeffs (Ynm*er)
shr = shAna(Fr); 

%Analyze u and v to obtain tangential coeffs 
snu = repmat(sin(u),1,d2); 

if nargin==1
    type='UV'; 
end

if strcmp(type,'GX')
    shu = shAna(Fu./snu); 
    shv = shAna(Fv./snu);
    % Build or load sparse matrix to convert to GYnm and Xnm coeffs
    T=LOCAL_get_VshMat(p); 
    Ctan = T*[shu;shv]; 

    shG = Ctan(1:sp,:); 
    shX = Ctan(sp+1:end,:);
    
    %Return
    shG(1)=0; shX(1)=0; 
    shG(p^2+1:sp,:)=0; shX(p^2+1:sp,:)=0; 
    Vshc=[shr;shG;shX]; 
elseif strcmp(type,'VW')
    shu = shAna(Fu./snu); 
    shv = shAna(Fv./snu);
    % Build or load sparse matrix to convert to GYnm and Xnm coeffs
    T=LOCAL_get_VshMat(p); 
    Ctan = T*[shu;shv]; 

    shG = Ctan(1:sp,:); 
    shX = Ctan(sp+1:end,:);
    
    % Convert CG and CX to CV and CW
    nn = floor(sqrt((1:sp)'-1)); 
    
    shV = repmat(-1./(2*nn+1),1,d2).*shr+repmat(nn./(2*nn+1),1,d2).*shG;
    shW = repmat(1./(2*nn+1),1,d2).*shr +repmat((nn+1)./(2*nn+1),1,d2).*shG;
    
    %Return V,W,X
    shW(1)=0; shX(1)=0; 
    shV(p^2+1:sp,:) = 0; shW(p^2+1:sp,:)=0; shX(p^2+1:sp,:)=0; 
    Vshc=[shV;shW;shX]; 
else
    % Standard: return r,u and v components
    shu = shAna(Fu.*snu); 
    shv = shAna(Fv.*snu); 
    Vshc=[shr;shu;shv]; 
end
    
end

function T = LOCAL_get_VshMat(p)

%TODO: add option to load

inm = @(n,m) n.^2+n+m+1; 
cv = @(m) 1i*m; 
cp1 = @(n,m) (n.*sqrt((n+1).^2-m.^2))./(sqrt((2*n+1).*(2*n+3)));
cm1 = @(n,m) -((n+1).*sqrt(n.^2-m.^2))./(sqrt((2*n+1).*(2*n-1)));
sp=(p+1)^2; 

ii = (1:sp)'; 
nn = floor(sqrt(ii-1)); 
mm=ii-nn.^2-nn-1;
dn = (nn.*(nn+1))+(nn==0); 
IGV = ii; JGV=ii+sp; VGV=-cv(mm)./dn; 
idm = nn>0 & nn>abs(mm);
IGU = ii(idm); JGU=ii(inm(nn(idm)-1,mm(idm))); 
VGU=cm1(nn(idm),mm(idm))./dn(idm);
IGU = [IGU;(1:p^2)']; 
JGU = [JGU;inm(nn(1:p^2)+1,mm(1:p^2))]; 
VGU=[VGU;cp1(nn(1:p^2),mm(1:p^2))./dn(1:p^2)];

IXU = ii+sp; JXU=ii; VXU=-VGV;
IXV=IGU+sp; JXV=JGU+sp; VXV=VGU; 

II=[IGU;IGV;IXU;IXV];
JJ=[JGU;JGV;JXU;JXV]; 
VV=[VGU;VGV;VXU;VXV];
T = sparse(II,JJ,VV,2*sp,2*sp); 

%TODO: add option to save

end

function testVShAna() 
p=6; 
sp=(p+1)^2; 
spm=p^2; 
np = 2*p*(p+1); 
n = 2; m=2; 
d2=1; 
inm = @(n,m) n.^2+n+m+1; 

[u,v]=gl_grid(p); 
Sc = SurfaceSph(shape_gallery(p,''));
% Normal Unit vector
er = [sin(u).*cos(v) sin(u).*sin(v) cos(u)];
% Tangential Unit vectors 
ev = [-sin(v) cos(v) zeros(np,1)]; 
eu = [cos(u).*cos(v) cos(u).*sin(v) -sin(u)];  

shrr = [rand(n^2,d2);zeros(sp-n^2,d2)]; fr = shSyn(shrr); 
shur = [rand(n^2,d2);zeros(sp-n^2,d2)]; fu = shSyn(shur); 
shvr = [rand(n^2,d2);zeros(sp-n^2,d2)]; fv = shSyn(shvr); 

snu = repmat(1./sin(u),1,d2); 

F = repmat(fr,3,1).*repmat(er(:),1,d2)+...
    repmat(snu.*fu,3,1).*repmat(eu(:),1,d2)+...
    repmat(snu.*fv,3,1).*repmat(ev(:),1,d2);

shuv = @(q) VshAna(q,'UV');
Vsh = shuv(F); %VshAna(F,'UV'); 
plot(log10(abs(Vsh)),'o'); hold on; plot(log10(abs([shrr;shur;shvr])),'or')

fprintf('\n Errors on shr, shu and shv\n')
display(norm(Vsh(1:sp,:)-shrr))
display(norm(Vsh(sp+1:2*sp,:)-shur))
display(norm(Vsh(2*sp+1:3*sp,:)-shvr))

%Test GX mode
shfgx = @(q) VshAna([q(1:3:end);q(2:3:end);q(3:3:end)],'GX');
YNr= Vnm('YNr',Sc,n,m,u,v); 
GY = Vnm('GY',Sc,n,m,u,v);
GY = reshape(reshape(GY,[np 3]).',[],1);
Xnm = Vnm('Xnm',Sc,n,m,u,v); 
F2 = YNr(:)+2*GY(:)+3*Xnm(:); 

Vsh2 = shfgx(F2); %VshAna(F2,'GX'); 
C=zeros(size(Vsh2)); C(inm(n,m))=1; C(inm(n,m)+sp)=2; C(inm(n,m)+2*sp)=3; 
fprintf('\n Errors on shr, shG and shX')
display(norm(Vsh2(1:spm,:)-C(1:spm,:)))
display(norm(Vsh2(sp+1:sp+spm,:)-C(sp+1:sp+spm,:)))
display(norm(Vsh2(2*sp+1:2*sp+spm,:)-C(2*sp+1:2*sp+spm,:)))

% Test VW mode
shfvw = @(q) VshAna([q(1:3:end);q(2:3:end);q(3:3:end)],'VW');
Vn= Vnm('Vnm',Sc,n,m,u,v); 
Wn = Vnm('Wnm',Sc,n,m,u,v); 

F3 = Vn(:)+2*Wn(:)+3*Xnm(:);
Vsh3 = shfvw(F3); % VshAna(F3,'VW'); 
fprintf('\n Errors on shV, shW and shX')
display(norm(Vsh3(1:spm,:)-C(1:spm,:)))
display(norm(Vsh3(sp+1:sp+spm,:)-C(sp+1:sp+spm,:)))
display(norm(Vsh3(2*sp+1:2*sp+spm,:)-C(2*sp+1:2*sp+spm,:)))

end
