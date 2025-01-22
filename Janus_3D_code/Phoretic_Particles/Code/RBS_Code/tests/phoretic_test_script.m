clear;clc;


fname = 'one_shell_workplz';
exp = 'fullsystem';
p = 8;

%C = [0 0.1 -1.5; 0 -0.1 1.5];
n = 6; dst=4; lx = 0:dst:dst*(n-1); lx=lx-mean(lx); 
[xx,yy,zz] = meshgrid(lx); 
%C = [xx(:) yy(:) zz(:)]; 
C = [-2 0 0];
C = C + 0.1*(2*rand(size(C))-1); 
n3 = size(C,1); 
%M0 = {eye(3), [1 0 0;0 1 0;0 0 -1]};
M0 = cell(n3,1); 
for i=1:n3
   r = rand(); 
   if r<0.5
      M0{i}=eye(3);  
   else
      M0{i}=[1 0 0;0 1 0;0 0 -1];  
   end
end

M0 = {-eye(3)};

mdist = 3;
tdisc = 'euler';
beta = 1.5;
flplot = 0;
% one particle speed is (4/45)*AC0*MC0
ct = sqrt(45/4); 
AC0 = ct; MC0=ct; 
Nt = 1000;
dt = 0.1;
A = @(theta) AC0*(1-cos(theta).^2);
M_theta = @(theta) MC0*cos(theta);

%set data for shell here
shell_data.p = 8;
shell_data.r = 10;
shell_data.Ldense = 0;
shell_data.Sdense = 0;



Test_Phoretic(fname,exp,p,C,M0,Nt,dt,mdist,tdisc,beta,flplot,A,M_theta, shell_data);
