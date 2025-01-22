function [Cx,Cy,Cz,Vx,Vy,Vz,Tx,Ty,Tz,tt] = RBS_extract_info(fname)

load(fname,'tt','Ct','VW'); 

N = size(Ct,1)-1;
ns = size(Ct{1},1); 

Cx = zeros(ns,1); Cy=Cx; Cz=Cx; Vx=Cx; Vy=Cx; Vz=Cx; Tx=Cx; Ty=Cx; Tz=Cx; 

for i=1:N
   if ~isempty(VW{i})
      Cx(:,i) = Ct{i}(:,1); 
      Cy(:,i) = Ct{i}(:,2);
      Cz(:,i) = Ct{i}(:,3);
      
      Vx(:,i) = VW{i}(1,:)'; 
      Vy(:,i) = VW{i}(2,:)';
      Vz(:,i) = VW{i}(3,:)';
      
      Tx(:,i) = VW{i}(4,:)';
      Ty(:,i) = VW{i}(5,:)';
      Tz(:,i) = VW{i}(6,:)';
   else
       N = i-1;
       break; 
   end
end

figure; 
plot3(Cx',Cy',Cz','-ob'); title('Particle centers'); 

display(N)
tt = tt(1:N); 
figure; plot(tt,Vx','-o'); title('Vx'); 
figure; plot(tt,Vy','-o'); title('Vy'); 
figure; plot(tt,Vz','-o'); title('Vz');

figure; plot(tt,Tx','-o'); title('Tx'); 
figure; plot(tt,Ty','-o'); title('Ty'); 
figure; plot(tt,Tz','-o'); title('Tz'); 


end