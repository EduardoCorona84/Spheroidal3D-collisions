video=true; 
shell = true;
vidname = 'one_controlled_flux';
clear X mn mx color; 
n3 = size(Ct{1},1); 
shell_p = 8;
nt = sum(tt>0);
np = 2*p*(p+1); 
warning off; 
idx = cell(n3,1);
color = cell(nt,n3); 


for i = 1:nt
    Concentration = Conc{i};
    shcol{i}=Concentration((shell_p+1)^2:end);
end


mxc = 0; 

for k=1:n3
   for i=1:nt
   p = (sqrt(1+2*(size(Xt{i},1)/n3))-1)/2; 
   np = 2*p*(p+1); 
   idx{k} = round(np*(k-1)+(1:np)); 
   if ~isempty(Xt{i})
      X{i,k}   = Xt{i}(idx{k},:)+repmat(Ct{i}(k,:),np,1); 
   if exist('phi')
      if ~isempty(phi{i}) 
        px = phi{i}(3*np*(k-1)+(1:3:3*np));
        py = phi{i}(3*np*(k-1)+(2:3:3*np));
        pz = phi{i}(3*np*(k-1)+(3:3:3*np));
        color{i,k} = max(abs([px py pz]'))'; 
        mxc = max(mxc,max(color{i,k})); 
      end
   else
      color{i,k}=0.5*ones(np,1); mxc=1;  
   end
   end
   end
end

if mxc>0
for k=1:n3
for i=1:nt
    if ~isempty(color{i,k})
        color{i,k} = color{i,k}/mxc; 
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
diam = 2*max(max(Xt{1})-min(Xt{1})); 
if n3>1
   axvec = [mn(1,:)-2.5*diam; mx(1,:)+2.5*diam]; axvec=axvec(:)'; 
else
   axvec = [Ct{1}(1)-2.5*diam; Ct{1}(1)+2.5*diam; Ct{1}(2)-2.5*diam; Ct{1}(2)+2.5*diam; ...
       Ct{1}(3)-2.5*diam; Ct{1}(3)+2.5*diam]; axvec=axvec(:)'; 
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




if isnumeric(shcol)
    if numel(shcol)==1
    shcolor = shcol*ones(size(XS,1),1);
    end
else
    shcolor = shcol{1}; 
end

alp = 0.25;   
vw = [-1 -1 0.5]; 
mxz = axvec(6); 
zoom = false;
shrd = 10;
if shell
   Sc = SurfaceSph(shape_gallery(shell_p,''));
   XS = shrd*reshape(Sc.cart.to_array,[],3);

   plotb(XS(:),real(shcolor)); 

   if ~zoom
   sh = shrd+0.5; 
   axvec = [-sh sh -sh sh -sh sh]; 
   axis(axvec); 
   end
   alpha(alp)
end

%caxis([0,1]); 
view(vw)
hold off; 

if video
set(gca,'nextplot','replacechildren'); 
set(gcf,'Renderer','zbuffer'); 
end

stride = 5;

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
   if ~isnumeric(shcol)
      shcolor = shcol{i};  
   end
   
   plotb(XS(:),real(shcolor)); 
   if ~zoom
   axvec = [-sh sh -sh sh -sh sh];
   end
   alpha(alp)
end

axis(axvec); %caxis([0,1]);
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