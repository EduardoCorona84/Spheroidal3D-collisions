
%video=true;
video=false; 
%phi=sigma;
zoom=true;
shell=false;
stride=10;
shrd=40;

%vidname='prince_test_control';
clear X mn mx color; 
n3 = size(Ct{1},1); 
p = (sqrt(1+2*(size(Xt{1},1)/n3))-1)/2; 
nt = length(Ct); 
np = 2*p*(p+1); 
warning off; 
idx = cell(n3,1);
color = cell(nt,n3); 

mxc = 0; 

for k=1:n3
   for i=1:nt
   p = (sqrt(1+2*(size(Xt{i},1)/n3))-1)/2; 
   np = 2*p*(p+1); 
   idx{k} = round(np*(k-1)+(1:np)); 
   if ~isempty(Xt{i})
      X{i,k}   = Xt{i}(idx{k},:)+repmat(Ct{i}(k,:),np,1); 
   if exist('phi')
      if ~isempty(phi{i}) && length(phi{i})==3*np*n3
        px = phi{i}(3*np*(k-1)+(1:3:3*np));
        py = phi{i}(3*np*(k-1)+(2:3:3*np));
        pz = phi{i}(3*np*(k-1)+(3:3:3*np));
        color{i,k} = max(abs([px py pz]'))';
        mxc = max(mxc,max(color{i,k}));
      elseif ~isempty(phi{i})
        muL = phi{i}(np*(k-1)+(1:np));  
        psiL = phi{i}(np*(n3+k-1)+(1:np));  
        color{i,k} = abs(muL); 
        mxc = max(mxc,max(color{i,k}));
      else
      color{i,k}=0.5*ones(np,1); %mxc=1;  
      end
   else
      color{i,k}=0.5*ones(np,1); %mxc=1; 
   end
   end
   end
end
display(mxc)
pause; 

maxc=1; 
if maxc>0
for k=1:n3
for i=1:nt
    if ~isempty(color{i,k})
        color{i,k} = color{i,k}/maxc; 
    end
end
end
end

nt = max(size(X,1)-1,1); 

for i=1:nt
mn(i,:) = min(Ct{i});
mx(i,:) = max(Ct{i});
end

%mn = min(mn); mx = max(mx); 
diam = max(max(Xt{1})-min(Xt{1})); 
if n3>1
   axvec = [mn(1,:)-1.5*diam; mx(1,:)+1.5*diam]; axvec=axvec(:)'; 
else
   axvec = [Ct{1}(1)-1.5*diam; Ct{1}(1)+1.5*diam; Ct{1}(2)-1.5*diam; Ct{1}(2)+1.5*diam; ...
       Ct{1}(3)-1.5*diam; Ct{1}(3)+1.5*diam]; axvec=axvec(:)'; 
end

if video
writerObj = VideoWriter([vidname '.avi']); 
writerObj.FrameRate=10; 
open(writerObj); 
end

figure; 
plotb(X{1,1}(:),color{1,1});
hold on;
for k=2:n3
plotb(X{1,k}(:),color{1,k});
end
axis(axvec);

alp = 0.3;  
shcol = 0.2; 
vw = [-1 0.5 0.5]; 
mxz = axvec(6); 

if shell
   Sc = SurfaceSph(shape_gallery(16,''));
   XS = shrd*reshape(Sc.cart.to_array,[],3);
   plotb(XS(:),shcol*ones(size(XS,1),1)); 
   if ~zoom
   sh = shrd+0.5; 
   axvec = [-sh sh -sh sh -sh mxz]; 
   axis(axvec); 
   end
   alpha(alp)
end

caxis([0,mxc]); 
display(mxc)
view(vw)
hold off; 

if video
set(gca,'nextplot','replacechildren'); 
set(gcf,'Renderer','zbuffer'); 
end


for i=1+stride:stride:nt-1
plotb(X{i,1}(:),color{i,1});
hold on;
for k=2:n3
plotb(X{i,k}(:),color{i,k});
end

if n3>1
axvec = [mn(i,:)-1.5*diam; mx(i,:)+1.5*diam]; axvec=axvec(:)'; 
else
axvec = [Ct{i}(1)-1.5*diam; Ct{i}(1)+1.5*diam; Ct{i}(2)-1.5*diam; Ct{i}(2)+1.5*diam; ...
       Ct{i}(3)-1.5*diam; Ct{i}(3)+1.5*diam]; axvec=axvec(:)';     
end

%q3(Xtr{i},Yv{i});
if shell
   plotb(XS(:),shcol*ones(size(XS,1),1)); 
   if ~zoom
   axvec = [-sh sh -sh sh -sh mxz];
   end
   alpha(alp)
end

axis(axvec); 
caxis([0,mxc]);
view(vw) 

if video
    
    frame=getframe(gcf); 
    writeVideo(writerObj,frame); 
else
   pause(0.1); 
end

hold off; 
end 

if video
close(writerObj);
end
