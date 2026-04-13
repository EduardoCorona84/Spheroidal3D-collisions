function F = interpvsh(F,pold,pnew)

%np = 2*pold*(pold+1); 
%nF = size(F,2); 
%F = reshape(interpsh(reshape(F,np,[]),pnew),[],nF); 

d2 = size(F,2);
np = @(p) 2*p*(p+1); 
sp = @(p) (p+1)^2; 

%{
d2 = size(F,2); 
np = @(p) 2*p*(p+1); 
sp = @(p) (p+1)^2; 

% Setup 
[u,v]=gl_grid(pold); 
% Normal Unit vector
er = [sin(u).*cos(v) sin(u).*sin(v) cos(u)];
% Tangential Unit vectors 
ev = [-sin(v) cos(v) zeros(np(pold),1)]; 
eu = [cos(u).*cos(v) cos(u).*sin(v) -sin(u)];  

dotp = @(x,y,np) x(1:np,:).*conj(y(1:np,:)) + x(np+1:2*np,:).*conj(y(np+1:2*np,:)) + x(2*np+1:3*np,:).*conj(y(2*np+1:3*np,:));

%Separate F into normal, u and v components
Fr = dotp(F,repmat(er(:),1,d2),np(pold)); 
Fv = dotp(F,repmat(ev(:),1,d2),np(pold)); 
Fu = dotp(F,repmat(eu(:),1,d2),np(pold));

%Apply interpsh
Frn = interpsh(Fr,pnew); 
Fvn = interpsh(Fv,pnew); 
Fun = interpsh(Fu,pnew);


% Setup 
[u,v]=gl_grid(pnew); 
% Normal Unit vector
er = [sin(u).*cos(v) sin(u).*sin(v) cos(u)];
% Tangential Unit vectors 
ev = [-sin(v) cos(v) zeros(np(pnew),1)]; 
eu = [cos(u).*cos(v) cos(u).*sin(v) -sin(u)];

%Form new F 
%F = Fr*er+Fu*eu+Fv*ev
F = repmat(Frn,3,1).*repmat(er(:),1,d2)+...
    repmat(Fun,3,1).*repmat(eu(:),1,d2)+...
    repmat(Fvn,3,1).*repmat(ev(:),1,d2);
%}

%{
[u,v]=gl_grid(pnew);

%Analyze r component to obtain radial coeffs (Ynm*er)
shr = shAna(Frn); 

%Analyze u and v to obtain tangential coeffs 
snu = repmat(sin(u),1,d2);
shu = shAna(Fun./snu); 
shv = shAna(Fvn./snu);

plot(log10(abs([shr;shu;shv])),'o')
F = VshSyn([shr;shu;shv],'UV'); 
    %}

F2 = [F(1:3:end,:);F(2:3:end,:);F(3:3:end,:)]; 
Fh = VshAna(F2,'VW'); 
zr = zeros(sp(pnew)-sp(pold),d2); 
Fhn = [Fh(1:sp(pold),:);zr;Fh(sp(pold)+1:2*sp(pold),:);zr;Fh(2*sp(pold)+1:3*sp(pold),:);zr];
F = VshSyn(Fhn,'VW'); 
F = reshape(reshape(F,[],3).',[],1);

end