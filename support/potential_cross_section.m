clear;clc;
load('tet_def');
lambda=1;
p=4;
np=2*p*(p+1);
type = 1;
%step=1;
Step=[1, sum(tt>0)];


for k=1:2


step = Step(k);
z = -2.055;
resolution = 0.01;
Centers=Ct{step};
r=ones(size(Centers,1),1);
Nodes=zeros(size(Centers,1)*np,3);
PL = psi_Lap{step}(1:end-2);
for i=1:size(Centers,1)
    range=np*(i-1)+1:np*i;
    for k=1:3
        Nodes(range,k)=Xt{step}(range,k)+Centers(i,k);
    end
end



xmin=min(Nodes(:,1)); xmax=max(Nodes(:,1));
ymin=min(Nodes(:,2)); ymax=max(Nodes(:,2));
[TrgX, TrgY]=meshgrid((xmin-1:resolution:xmax+1),(ymin-1:resolution:ymax+1));
Targ=[TrgX(:),TrgY(:),z*ones(size(TrgX(:)))];



SLparams=Setparams(p,Centers,r,'SL_LMOD_3D','',1,1,5,1,lambda); SLparams.X=Nodes;
DLparams=Setparams(p,Centers,r,'DL_LMOD_3D','',1,1,5,1,lambda); DLparams.X=Nodes;
SL = @(v,T) VSh_Mod_MatVec_RB_trg(v,T,T,SLparams);
DL = @(v,T) VSh_Mod_MatVec_RB_trg(v,T,T,DLparams);


if (type == 1) %Amphiphilic
    mu = real(PL(1:end/2));
    phi = real(SL(mu,Targ)+ DL(mu,Targ));
    
elseif (type == 2)
    mu = real(PL{step}(1:end/2));
    psi = real(psi_Lap{step}(end/2+1:end));
    E=[-1 0 0];
    phi = real(SL(psi,Targ) + DL(mu,Targ));%  + Targ*E';
    
end

figure;
c = jet(256);
phi = reshape(phi,size(TrgX,1),[]);

dphi=zeros(size(phi));
for i=2:size(phi,1)-1
   dphi(i,:) = (phi(i-1,:) - phi(i+1,:))/(2*resolution); 
    
end
dphi((abs(dphi) > 1)) = 0;


surf(TrgX, TrgY, phi ,'EdgeColor', 'None', 'facecolor', 'interp');
%surf(TrgX, TrgY, dphi ,'EdgeColor', 'None', 'facecolor', 'interp');

colormap(c);
 cb =  colorbar;
 cb.LineWidth = 1.5;

shading interp;
view(2);
axis equal; 
axis off;
%{
[dx dy]=gradient(-phi,resolution);
hold on
%quiver(TrgX,TrgY,dx,dy);
[sTrgX, sTrgY]=meshgrid((xmin-1:10*resolution:xmax+1),(ymin-1:10*resolution:ymax+1));
streamline(TrgX,TrgY,dx,dy,sTrgX,sTrgY);
hold off
%}


end