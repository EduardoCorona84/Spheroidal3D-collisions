clear;
p=24;
lambda=1e-10;
rng('default');
type='DMat';


sp = @(p) (p+1)^2;  %(number of harmonic coefficients);
np = @(p) 2*(p+1)*p; %(number of quadrature points);
ii=(1:sp(p))';
nn=floor(sqrt(ii-1)); %indices of n
%generate random seeds


swgt=rand(sp(p),1);

for i=1:sp(p)        %exponentially decaying weights
    swgt(i)=swgt(i)*(1/2)^(floor(sqrt(i-1)));
end
swgt(sp(p-1)+1:end)=0;


%{
swgt=zeros(sp(p),1);
swgt(2)=1;
%}
  
%}
%evaluates via Sh_Kernel_Eval
Sc = SurfaceSph(shape_gallery(p,''));
X =reshape(Sc.cart.to_array,[],3); 

Ctrg = [0 0 5;5 0 0;0 5 0;-5 0 0;0 -5 0;0 0 -5]; 
Rtrg=2; 

for i=1:size(Ctrg,1)
Trg= SurfaceSph(shape_gallery(p,''));
if i==1
    XT =Rtrg*reshape(Trg.cart.to_array,[],3)+repmat(Ctrg(1,:),np(p),1); 
else
    XT=vertcat(XT,Rtrg*reshape(Trg.cart.to_array,[],3)+repmat(Ctrg(i,:),np(p),1));
end
end
% v = th in [0,2*pi], u = phi in [0,pi] with 0 -> north pole
[v,u, rad]=cart2sph(XT(:,1),XT(:,2),XT(:,3));
v(v<0)=v(v<0)+2*pi; u=pi/2-u;
%}

K=Numerical_Kernel(XT,Sc,p,lambda,type);
Y=shSyn(swgt);
KY=K*Y;


rad=1;
KS=(Sh_Kernel_Eval_off(swgt,type,0,1,rad,[],[],[]));
tic
KSM=((Sh_Mod_Kernel_Eval_off(swgt,type,0,1,rad,[],[],[],lambda)));

Error(1)=log10(norm((KS-KSM))/norm(KS));
%Error(2)=log10(norm(KSM-KY)/norm(KY));
fprintf('\n-------------------------------------------------------------'); 
fprintf('\n---------------Modified Laplace lambda=%2.2g-----------------',lambda);
fprintf('\n---------------Rtrg = %2.2g, Ctrg=[%1.1g,%1.1g,%1.1g]-----------------------',Rtrg,Ctrg(1),Ctrg(2),Ctrg(3));
fprintf('\n Rel error ||KLap - Kmod||/||KLap|| = %4.4g',Error(1));
%fprintf('\n Rel error ||Kmodnum - Kmod||/||Kmodnum|| = %4.4g',Error(2)); 
fprintf('\n-------------------------------------------------------------\n');
%display(Error);




