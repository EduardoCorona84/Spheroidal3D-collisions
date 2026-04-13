
alpha=0.02;
rng('default');
n3=size(Ct{1},1);
xx=Ct{1}(:,1);
timestep = sum(tt>0);
%v=[ -2*alpha*xx(:), zeros(size(xx(:))), ones(size(xx(:))) ];
%init_dir=[ v; -v];
%init_dir=zeros(n3,3);
%init_dir = repmat([0 1 0],n3,1);

init_dir =   [  0.6294   -0.2155   -0.7620
    0.8116    0.3110   -0.0033
   -0.7460   -0.6576    0.9195
    0.8268    0.4121   -0.3192
    0.2647   -0.9363    0.1705
   -0.8049   -0.4462   -0.5524
   -0.4430   -0.9077    0.5025
    0.0938   -0.8057   -0.4898
    0.9150    0.6469    0.0119
    0.9298    0.3897    0.3982
   -0.6848   -0.3658    0.7818
    0.9412    0.9004    0.9186
    0.9143   -0.9311    0.0944
   -0.0292   -0.1225   -0.7228
    0.6006   -0.2369   -0.7014
   -0.7162    0.5310   -0.4850
   -0.1565    0.5904    0.6814
    0.8315   -0.6263   -0.4914
    0.5844   -0.0205    0.6286
    0.9190   -0.1088   -0.5130
    0.3115    0.2926    0.8585
   -0.9286    0.4187   -0.3000
    0.6983    0.5094   -0.6068
    0.8680   -0.4479   -0.4978
    0.3575    0.3594    0.2321
    0.5155    0.3102   -0.0534
    0.4863   -0.6748   -0.2967];

%{
for i=1:n3
init_dir(i,:) = (Mt{timestep,i}'*[-1 0 0]')';
end
%}
%init_dir = repmat([1 0 0],n3,1);
%init_dir=2*(rand(n3,3)-0.5*ones(n3,3));
%init_dir=[1 0 0;-1 0 0];
%init_dir=repmat([0 -1 0],n3,1);
%init_dir=[repmat([0 0 1],n3/2,1); repmat([0,0,-1],n3/2,1)];
%init_dir=[repmat(v,n3/2,1); repmat(-v,n3/2,1)]; init_dir=init_dir(1:n3,:);
%init_dir=-Ct{1};
init_dir=init_dir./repmat(sqrt(init_dir(:,1).^2+init_dir(:,2).^2+init_dir(:,3).^2),1,3);


file='prince_test_FMM_lattice_deux';

Current_Dir=zeros(size(init_dir));

for j=1:5:timestep
    Cinit=Ct{j};
    for i=1:n3
      
       %Center=Fparams.parmod.C(i,:);
       %Xsphere=Xsphere-repmat(Center,np,1);
        
       
       filename=strcat(file,int2str((j-1)/5),'.vtk');
        Current_Dir(i,:)=Mt{j,i}*init_dir(i,:)';
       %fprintf('\n max surf label = %e, min surf label = %e ',max(surfacelabel(:,i)),min(surfacelabel(:,i)));
       %plot(surfacelabel(:,i)); 
       
        vtkwrite(filename, 'structured_grid',Cinit(:,1), Cinit(:,2), Cinit(:,3),'vectors','dir',Current_Dir(:,1),Current_Dir(:,2),Current_Dir(:,3));
    end
end 
    
    
    
    
    
    