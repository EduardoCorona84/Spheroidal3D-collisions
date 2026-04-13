function F = VshSyn(Vshc,type,isReal)
%{
F = VSHSYN(SHC,type,isReal) - Calculates the inverse vector spherical 
transform for the given set of vector spherical harmonics coefficients SHC. 
Each column of SHC is the coefficients of a vector field. 
The kind of vector spherical harmonics used is indicated by type:

type='UV' -> Vshc =[shr;shu;shv] analysis in radial, u and v directions
type='GX' -> Vshc =[shr;shG;shX] basis functions {Ynm*Nr,GYnm,Xnm}
type='VW' -> Vshc =[shV;shW;shX] basis functions {Vnm,Wnm,Xnm}
 
isReal is the flag that indicates
that the desired function F is real valued.

SEE ALSO: VSHANA, VNM, YNM.
%}

%-- Checking the input
if(nargin==0), testVShSyn(); return; end
if(nargin<3), isReal = false; end
if(nargin<2), type = 'UV'; end

% Sizes and parameters
[d1,d2] = size(Vshc);
sp = d1/3; %sp = (p+1)^2 # of basis elements
p = sqrt(sp)-1; %d1 = 3(p+1)^2
np = 2*p*(p+1); 
[u,v]=gl_grid(p); 
% Normal Unit vector
er = [sin(u).*cos(v) sin(u).*sin(v) cos(u)];
% Tangential Unit vectors 
ev = [-sin(v) cos(v) zeros(np,1)]; 
eu = [cos(u).*cos(v) cos(u).*sin(v) -sin(u)];  

% Recover F in spherical coordinates (Fr,Fu,Fv)
if strcmp(type,'UV')
    %Apply scalar synthesis to each component of Vshc
    Fcomp = shSyn([Vshc(1:sp,:) Vshc(sp+1:2*sp,:) Vshc(2*sp+1:3*sp,:)],isReal); 
elseif strcmp(type,'GX')
    % Fcomp = [Fr,Fg,Fx]
    shr = Vshc(1:sp,:); 
    shG = Vshc(sp+1:2*sp,:); 
    shX = Vshc(2*sp+1:3*sp,:);
    
    % [shu;shv] = T*[shG;shX]
    T = LOCAL_get_VshMat(p);
    shuv = T*[shG;shX]; 
    
    Fcomp = shSyn([shr shuv(1:sp,:) shuv(sp+1:2*sp,:)],isReal); 
elseif strcmp(type,'VW')
    % Fcomp = [cV,cW,cX]
    shV = Vshc(1:sp,:); 
    shW = Vshc(sp+1:2*sp,:); 
    shX = Vshc(2*sp+1:3*sp,:);
    
    % Get Fr and Fg from Fv and Fw
    nn = floor(sqrt((1:sp)'-1)); 
    % V = G-(n+1)YNr, W=G+nYNr
    shr = repmat(-(nn+1),1,d2).*shV+repmat(nn,1,d2).*shW;
    shG = shV+shW;
    
    % [shu;shv] = T*[shG;shX]
    T = LOCAL_get_VshMat(p);
    shuv = T*[shG;shX]; 
    
    Fcomp = shSyn([shr shuv(1:sp,:) shuv(sp+1:2*sp,:)],isReal); 
else
    fprintf('\n Type not recognized \n')
    F=[]; return; 
end

% Fcomp = [Fr,Fu/sin(u),Fv/sin(u)]
snu = repmat(sin(u),1,d2); 
Fr = Fcomp(:,1:d2); 
Fu = Fcomp(:,d2+1:2*d2)./snu; 
Fv = Fcomp(:,2*d2+1:3*d2)./snu; 

%F = Fr*er+Fu*eu+Fv*ev
F = repmat(Fr,3,1).*repmat(er(:),1,d2)+...
    repmat(Fu,3,1).*repmat(eu(:),1,d2)+...
    repmat(Fv,3,1).*repmat(ev(:),1,d2);

end

function T = LOCAL_get_VshMat(p)

%TODO: add option to load

%Form T matrix if not available
inm = @(n,m) n.^2+n+m+1; 
cv = @(m) 1i*m; 
cp1 = @(n,m) (n.*sqrt((n+1).^2-m.^2))./(sqrt((2*n+1).*(2*n+3)));
cm1 = @(n,m) -((n+1).*sqrt(n.^2-m.^2))./(sqrt((2*n+1).*(2*n-1)));
sp=(p+1)^2; 

ii = (1:sp)'; 
nn = floor(sqrt(ii-1)); 
mm=ii-nn.^2-nn-1;
IGV = ii; JGV=ii+sp; VGV=cv(mm); 
idm = nn>0 & nn>abs(mm);
IGU = ii(idm); JGU=ii(inm(nn(idm)-1,mm(idm))); 
VGU=cm1(nn(idm),mm(idm));
IGU = [IGU;(1:p^2)']; 
JGU = [JGU;inm(nn(1:p^2)+1,mm(1:p^2))]; 
VGU=[VGU;cp1(nn(1:p^2),mm(1:p^2))];

IXU = ii+sp; JXU=ii; VXU=-VGV;
IXV=IGU+sp; JXV=JGU+sp; VXV=VGU; 

II=[IGU;IGV;IXU;IXV];
JJ=[JGU;JGV;JXU;JXV]; 
VV=[VGU;VGV;VXU;VXV];
T = sparse(JJ,II,VV,2*sp,2*sp); 

end

function [T,Jrow,Jcol] = LOCAL_get_InvVshMat(p)

%TODO: add option to load

%(1) Form T matrix if not available
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

%(2) Eliminate rows and columns and find good ordering to invert
T([1 sp+1],:)=[]; 
npm0=p^2+p+1; 
T(:,[npm0 npm0+sp])=[];

mmr=mm(2:end); mmc=mm; mmc(npm0)=[];
Jr = find(mmr==0); Jc = find(mmc==0); 
for i=1:p
Jr=[Jr; find(mmr==-i); find(mmr==i)]; 
Jc=[Jc; find(mmc==-i); find(mmc==i)];
end
Jrow = zeros(2*sp-2,1); Jcol=Jrow; 
Jrow(1:2:end) = Jr; Jrow(2:2:end)=Jr+sp-1;
Jcol(1:2:end) = Jc; Jcol(2:2:end)=Jc+sp-1;

%TODO: add option to save

end

function shuv=LOCAL_solveT(T,C,Jrow,Jcol)
d2=size(C,2); 
sp = size(C,1)/2; 
p = sqrt(sp)-1;  
npm0=p^2+p+1; 
shuv = zeros(2*sp,d2); 
idx = 1:2*sp; 
idi=(idx~=1 & idx~=sp+1)'; 
idj=(idx~=npm0 & idx~=npm0+sp)';

B = C(idi,:); 
tmp=zeros(2*sp-2,d2); 
tmp(Jcol,:) = T(Jrow,Jcol)\B(Jrow,:); 
shuv(idj,:)=tmp; 

end

function testVShSyn()

p=16; 
sp=(p+1)^2; 
np = 2*p*(p+1); 
n = 8; 
d2=1; 

[u,v]=gl_grid(p); 
% Normal Unit vector
er = [sin(u).*cos(v) sin(u).*sin(v) cos(u)];
% Tangential Unit vectors 
ev = [-sin(v) cos(v) zeros(np,1)]; 
eu = [cos(u).*cos(v) cos(u).*sin(v) -sin(u)];  

shrr = [rand(n^2,d2);zeros(sp-n^2,d2)]; fr = shSyn(shrr); 
shur = [rand(n^2,d2);zeros(sp-n^2,d2)]; fu = shSyn(shur); 
shvr = [rand(n^2,d2);zeros(sp-n^2,d2)]; fv = shSyn(shvr); 

%shrr = zeros(sp,1); shur=zeros(sp,1); shvr=zeros(sp,1); shur(10)=1; 
Vshrd = [shrr;shur;shvr]; 

snu = repmat(1./sin(u),1,d2); 

F = repmat(fr,3,1).*repmat(er(:),1,d2)+...
    repmat(snu.*fu,3,1).*repmat(eu(:),1,d2)+...
    repmat(snu.*fv,3,1).*repmat(ev(:),1,d2);

Vsh = VshAna(F,'UV'); 
figure; plot(real(F))
F2 = VshSyn(Vsh,'UV'); 

fprintf('\n Error reconstructing F from Vsh (UV)\n')
display(norm(F-F2))

Vshgx = Vshrd; %Vshgx(sp+1:2*sp,:) = 0; 
FG = LOCAL_VshSyn_brute(Vshgx,'GX'); 
FG2 = zeros(size(FG)); 
for i=1:d2
FG2(:,i) = [FG(1:3:end,i) ; FG2(2:3:end,i) ; FG(3:3:end,i)];
end
figure; plot(real(FG2),'r')
Vsh2 = VshAna(FG2,'GX'); 
FG3 = VshSyn(Vsh2,'GX'); 

fprintf('\n Error reconstructing F from Vsh (GX)\n')
display(norm(FG3-FG2))
figure; plot(log10(abs(Vshgx-Vsh2)),'or')

FV = LOCAL_VshSyn_brute(Vshrd,'VW'); 
FV2 = [FV(1:3:end) ; FV(2:3:end) ; FV(3:3:end)]; 
figure; plot(real(FV),'k')
Vsh3 = VshAna(FV2,'VW'); 
FV3 = VshSyn(Vsh3,'VW');
figure; plot(log10(abs(Vshrd-Vsh3)),'k')

fprintf('\n Error reconstructing F from Vsh (VW)\n')
display(norm(FV3-FV2))

end

function F = LOCAL_VshSyn_brute(Vsh,type)

[d1,d2] = size(Vsh);
sp = d1/3; 
p = sqrt(sp)-1;
np = 2*p*(p+1); 
[u,v] = gl_grid(p);
Sc = SurfaceSph(shape_gallery(p,''));
Nr = reshape(Sc.geoProp.nor.to_array,[],3);
F = zeros(3*np,1); 

if strcmp(type,'VW')
   V1 = 'Vnm'; V2 = 'Wnm'; V3 = 'Xnm';  
elseif strcmp(type,'GX')
   V1 = 'YNr'; V2 = 'GY'; V3 = 'Xnm';  
end

for k=0:p
   ind = k^2+1:(k+1)^2; 
   Vk = Vnm(V1,Sc,k,[],u,v,Nr); 
   Wk = Vnm(V2,Sc,k,[],u,v,Nr);
   
   if strcmp(V2,'GY')
      szn = 2*k+1; 
      Wk = reshape(permute(reshape(Wk,[np 3 szn]),[2 1 3]),[],szn);
   end
   
   Xk = Vnm(V3,Sc,k,[],u,v,Nr); 
   
   F = F + [Vk Wk Xk]*Vsh([ind' ; ind'+sp ; ind'+2*sp]); 
end

end