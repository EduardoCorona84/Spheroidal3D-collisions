% Curvature Stokian flow over one/large system of prolates/mixes

clear; close all;
ns_inside=30; p=4; np=2*p*(p+1);

%------ Set up system of spheroids -------------
if ns_inside==3
    ns=ns_inside+1;
    u0=[1.1 1.2 1.3 3/sqrt(7)];
    % ns=ns_inside;
    % u0=[1.1 1.2 1.3];
    pars=SpheroidalParameters;
    pars.matvec_eta = 10;
    pars.isReal=0;
    pars.u0=u0;
    pars.p=p;
    alist = 1./u0;
    alist(end) = 8./sqrt(u0(end)^2+1);
    pars.a=alist;
    pars.oblate = [0 0 0 1];
    % pars.oblate=[0 0 0];

    pars.centers = [0 0 0; 4 0 0; 1.6 1.6 1.6; 0 0 0];
    thetas = [0 pi/10 5*pi/3 0];
    phis = [0 0 pi/5 0];

    % pars.centers = [0 0 0; 4 0 0; 1.6 1.6 1.6];
    % thetas = [0 pi/10 5*pi/3];
    % phis = [0 0 pi/5];

    Ri=zeros(3,3,ns);
    for i=1:ns
        thetai=thetas(i);
        phii=phis(i);
        Riy=[cos(thetai) 0 sin(thetai); 0 1 0; -sin(thetai) 0 cos(thetai)];
        Riz=[cos(phii) -sin(phii) 0; sin(phii) cos(phii) 0; 0 0 1];
        Ri(:,:,i)=Riz*Riy;
    end
    pars.Rmat=Ri;

elseif ns_inside==1
    ns=ns_inside;
    u0=1.1;
    pars=SpheroidalParameters;
    pars.matvec_eta = 10;
    pars.isReal=0;
    pars.u0=u0;
    pars.p=p;
    pars.a=1/u0;
    pars.oblate=0;
    pars.centers=[0 0 0];
    thetas=0;
    phis=0;

    pars.Rmat=eye(3,3);

else
    [u0,a,oblate,centers,Rmats,ns_actual,dbox]=setup_large(ns_inside,p);
    ns=ns_actual+1;
    pars=SpheroidalParameters;
    pars.matvec_eta = 10;
    pars.isReal=0;
    pars.u0=[u0 3/sqrt(7)]; 
    pars.a = [a dbox*sqrt(7)/4];
    pars.oblate=[oblate 1];
    pars.centers=[centers; 0 0 0]; 
    pars.Rmat=cat(3,Rmats,eye(3));
    % pars.plot();
    % Xsph=pars.get_X();
    % for ii=1:ns_actual
    %     cur_sph=[Xsph((ii-1)*np+1:ii*np,1);Xsph((ii-1)*np+1:ii*np,2);Xsph((ii-1)*np+1:ii*np,3)];
    %     plotb(cur_sph); hold on;
    % end
    % hold off
    % TODO: or use sigma to shade? plotb(pars.get_X(),sigma);
end

%------------------------------------------------------
% ------- Writing setup to file -----------------------------------
file_name='./data/stokes_'+string(ns_inside)+'_ptcls_setup.txt';
fID=fopen(file_name,'w');
fprintf(fID,'p=%d\n',p); % changed: store p value as well.
fprintf(fID,"ns total=%d\n",ns); % changed: now stores total number of spheroids (including outer shell)

fprintf(fID,"centers:x,y,z\n");
fprintf(fID,'%f, %f, %f\r\n',pars.centers');
fprintf(fID,"\n ");

fprintf(fID,'u0s:\n');
fprintf(fID,'%f\r\n',pars.u0);
fprintf(fID,"\n ");

fprintf(fID,"a's:\n");
fprintf(fID,'%f\r\n',pars.a);
fprintf(fID,"\n ");

fprintf(fID,'oblate:\n');
fprintf(fID,'%i\r\n',pars.oblate);
fprintf(fID,"\n ");

fprintf(fID,'Rmat: x,y,z\n');
for i=1:ns
    Mat_temp = pars.Rmat(:,:,i);
    fprintf(fID,'%f, %f, %f\r\n',Mat_temp');
    fprintf(fID,'\n');
end
fclose(fID);
%----------------------------------------------------

% ------------- For curvature flow, sigma=K*n ---------------------
% %{
% figure;
tic
% [Y,~]=pars.get_X();
NrY=pars.get_Norm_rot(p);
sigma_x=zeros(np,1,ns); sigma_y=zeros(np,1,ns); sigma_z=zeros(np,1,ns);
pars.sigma=sigma_x;
pars.get_shc();
[Y,~]=pars.get_X();
[theta,phi]=gl_grid(p);
v=cos(theta);
for i=1:ns
    Yi=Y((i-1)*np+1:i*np,:);
    Nri=NrY((i-1)*np+1:i*np,:);
    ai=pars.a(i); 
    ui=pars.u0(i);
    if ~pars.oblate(i)
        K=ui.*(2.*ui.^2-v.^2-1)./(2*ai)./sqrt(ui.^2-1)./(ui.^2-v.^2).^(3/2);
    else
        K=ui.*(2.*ui.^2+v.^2+1)./(2*ai)./sqrt(ui.^2+1)./(ui.^2+v.^2).^(3/2);
    end
    sigma_x(:,:,i)=-K.*Nri(:,1);
    sigma_y(:,:,i)=-K.*Nri(:,2);
    sigma_z(:,:,i)=-K.*Nri(:,3);
    sigma_norm=sqrt((K.*Nri(:,1)).^2+(K.*Nri(:,2)).^2+(K.*Nri(:,3)).^2);
    if i~=ns
        plotb([Yi(:,1);Yi(:,2);Yi(:,3)],sigma_norm); hold on;
    else
        scatter3(Yi(:,1),Yi(:,2),Yi(:,3),20.*ones(size(Yi(:,1))),sigma_norm,'filled'); hold on;
    end
end
tend=toc;
fprintf("\n time to calculate curvature as density: %d seconds.",tend);

title("Curvature norm")
colorbar;
hold off;

pic_name='./data/stokes_'+string(ns_inside)+'_ptcls_sigmaNorm.png';
fID=fopen(pic_name,'w');
saveas(gcf,pic_name);
fclose(fID);
% %}
%---------------------------------------------

%----- Set up domain points for velocity field --------------------
% %{
tic
DX=(-dbox:0.51:dbox+.08)'; DY=DX; DZ=DX;
[XX,YY,ZZ]=meshgrid(DX,DY,DZ);
Xbox=[XX(:),YY(:),ZZ(:)];

is_clear_bound=check_clear(pars.centers(end,:),pars.u0(end),pars.a(end),pars.oblate(end),Xbox);
Xgrid=Xbox(is_clear_bound==0,:);

% only look at points that are also outside the interior spheroids.
clear_mat=ones(size(Xgrid,1),1);
for i=1:ns-1 % check overlap with every particle inside outer constraint.
    is_clear=check_clear(pars.centers(i,:),pars.u0(i),pars.a(i),pars.oblate(i),Xgrid);
    clear_mat=clear_mat.*is_clear; % if 0 (not clear) for one spheroid i, will be 0 at the end.
end
Xeval=Xgrid(clear_mat==1,:);

Stk_x=zeros(size(Xgrid,1),1); Stk_y=Stk_x; Stk_z=Stk_x;
tend=toc;
fprintf("\n time to set up domain points: %d seconds.",tend);
% %}
%--------------------------------------------------------------------------

%--------- Stokes Layer potential Evaluation ------------------------------
% %{
tic;
[Stk_x_eval,Stk_y_eval,Stk_z_eval]=L2Stk(Xeval,pars,sigma_x,sigma_y,sigma_z,ns,theta,phi);
tend=toc;
fprintf("\n time to perform Stokes calculation for p=%d: %d seconds.", p, tend);

Stk_x(clear_mat==1,:)=real(Stk_x_eval);
Stk_y(clear_mat==1,:)=real(Stk_y_eval);
Stk_z(clear_mat==1,:)=real(Stk_z_eval);

Stk_XX=zeros(size(XX)); Stk_YY=Stk_XX; Stk_ZZ=Stk_XX;
[nx,ny,nz]=size(XX);
stklist_ind=1;
for zz=1:nz
    for yy=1:ny
        for xx=1:nx
            box_ind=ny*nx*(zz-1)+nx*(yy-1)+xx;
            if is_clear_bound(box_ind)==0 % if grid point inside outer constraint
                Stk_XX(xx,yy,zz)=Stk_x(stklist_ind);
                Stk_YY(xx,yy,zz)=Stk_y(stklist_ind);
                Stk_ZZ(xx,yy,zz)=Stk_z(stklist_ind);
                stklist_ind=stklist_ind+1;
            end
        end
    end
end
% %}
%--------------------------------------------------------------------------

% ------ Write Stokes matrices to file. -------------------------------------
% %{
file_name='./data/stokes_'+string(ns_inside)+'_ptcls_x.txt';
fID=fopen(file_name,'w');
fprintf(fID,'%f %f %f\r\n',Stk_XX);
fclose(fID);

file_name='./data/stokes_'+string(ns_inside)+'_ptcls_y.txt';
fID=fopen(file_name,'w');
fprintf(fID,'%f %f %f\r\n',Stk_YY);
fclose(fID);

file_name='./data/stokes_'+string(ns_inside)+'_ptcls_z.txt';
fID=fopen(file_name,'w');
fprintf(fID,'%f %f %f\r\n',Stk_ZZ);
fclose(fID);
% %}

%----2D streamline plotting------------------------------------------------
% %{
figure;
steps=floor(length(DZ)/12);
slices=(length(DZ)-12*steps)/2+steps.*(1:12);
% slices=[10,12,14,16,18,20,22,24];
nslice=length(slices);
for s=1:nslice
    subplot(3,4,s)
    ZZ_slice=ZZ(:,:,slices(s));
    title("Flow line for z="+string(ZZ_slice(1,1)));
    HS = streamline(stream2(XX(:,:,slices(s)), YY(:,:,slices(s)), Stk_XX(:,:,slices(s)), Stk_YY(:,:,slices(s)), XX(:,:,slices(s)), YY(:,:,slices(s))));
    % FFx, FFy are the velocities at XX, YY points (format should be as in meshgrid)
    % set(HS, 'LineWidth', 0.02, 'Color', 'cyan'); axis off;  hold on
    axis equal;
    set(HS, 'LineWidth',1); 
    xlabel("x");ylabel("y");
    % hold on; plot([Xbox(:,1); Xbox(1,1)], [Xbox(:,2); Xbox(1,2)],'r');
end
hold off;
pic_name='./data/stokes_'+string(ns_inside)+'_ptcls_streamline.png';
fID=fopen(pic_name,'w');
saveas(gcf,pic_name);
fclose(fID);
%--------------------------------------------------------------------------

% %---------2D quiver plot projection-------------
% figure;
% slices=[4,6,10,14,18,22,26,28];
% nslice=length(slices);
% for s=1:nslice
%     subplot(2,4,s)
%     ZZ_slice=ZZ(:,:,slices(s));
%     title("Velocity field for z="+string(ZZ_slice(1,1)));
%     HS = quiver(XX(:,:,slices(s)), YY(:,:,slices(s)), Stk_XX(:,:,slices(s)), Stk_YY(:,:,slices(s)));
%     axis equal;
%     set(HS, 'LineWidth',1); 
%     xlabel("x");ylabel("y");
% end
% hold off;
% pic_name='./data/stokes_'+string(ns_inside)+'_ptcls_field.png';
% fID=fopen(pic_name,'w');
% saveas(gcf,pic_name);
% fclose(fID);
% % %}
%--------------------------------------------------------------------------


%% FUNCTIONS
function is_clear = check_clear(cntr,u0,a,ob,xtrg)
    % check if points xtrg are inside a given spheroid

    % input:
    %       cntr: 1 x 3 position of spheroid.
    %       u0, a, ob: parameters of the spheroid i.
    %       xtrg: n x 3 points to determine.

    % output: 
    %       is_clear: n x 1 boolean of whether each point is "clear", i.e.
    %       outside the given spheroid.

    xdis=xtrg(:,1)-cntr(1);
    ydis=xtrg(:,2)-cntr(2);
    zdis=xtrg(:,3)-cntr(3);
    if ob
        formula=(xdis.^2+ydis.^2)./a^2./(u0^2+1)+zdis.^2./a^2./u0^2;
    else
        formula=(xdis.^2+ydis.^2)./a^2./(u0^2-1)+zdis.^2./a^2./u0^2;
    end
    is_clear=formula>1;
end

function [u0,a,oblate,centers,Rmats,ns_actual,dbox]=setup_large(ns,p)
    succ_flag=0;
    dbox=8;

    % to avoid long wait time on MaximumClique, catch time out messages,
    % then regenerate a system of spheroids using a larger otuer shell.
    % If still timed out, increase the major axis of the outer shell by 9/8, 
    % until Max Clique was achieved without going over time.

    % TODO later? Instead of increasing outer shell, decrease radii of
    % inside spheres at random.

    while ~succ_flag
        try
            [sys_params,inds_out,ns_actual] = call_randspheroids(ns,p,dbox);
            succ_flag=1;
        catch ME
            fprintf("\n timed out, retry generation with larger box.");
            dbox=dbox*9/8;
        end
    end

    fprintf("\n wanted a system of %d particles, resulted in %d particles.\n", ns, ns_actual);

    centers=zeros(ns_actual,3);
    Rmats=zeros(3,3,ns_actual);
    oblate=zeros(1,ns_actual);
    u0=zeros(1,ns_actual);
    a=zeros(1,ns_actual);

    for ii2=1:ns_actual
        ind2=inds_out(ii2);
        centers(ii2,:)=sys_params(ind2).C;
        Rmats(:,:,ii2)=sys_params(ind2).R;
        % % from when only two shapes 'ellipseZ' or 'oblateZ'.
        % oblate(ii2)=strcmp(sys_params(ii2).shape,'oblateZ');

        Aii=sys_params(ind2).a;
        Cii=sys_params(ind2).c;
        if Aii>Cii
            % R = major/minor, Aii/Cii for oblates.
            oblate(ii2)=1;
            u0(ii2) = 1/sqrt((Aii/Cii)^2-1);
            a(ii2) = 1/sqrt(u0(ii2)^2+1);
        else
            % R = major/minor, Cii/Aii for prolates.
            oblate(ii2)=0;
            u0(ii2) = 1/sqrt(1-(Aii/Cii)^2);
            a(ii2) = 1/u0(ii2);
        end
    end

    % % from when only two shapes 'ellipseZ' or 'oblateZ'.
    % u0=2/sqrt(3).*ones(size(oblate));
    % a=~oblate./u0 + oblate./sqrt(u0.^2+1);

end

function pcp = PtChargePotential(ptch,Xptch,Y)
M=length(ptch);
% np=length(Y);
np=size(Y,1);
Rptch=zeros(np,M);
for i=1:M
    Rptch(:,i)=sqrt(sum((Xptch(i,:)-Y).^2,2));
end
pcp=1./(4*pi*Rptch)*ptch; 
end

function flux=PtChargeFlux(ptch,Xptch,Y,NrY)
    M=length(ptch);
    np=length(Y);
    Rptch=zeros(np,M);
    EdotN=zeros(np,M);
    for i=1:M
        r=Xptch(i,:)-Y;
        Rptch(:,i)=sqrt(sum(r.^2,2));
        EdotN(:,i)=dot(NrY,r./Rptch(:,i),2);
    end
    flux=EdotN./(4*pi*Rptch.^2)*ptch;
end

function E=PtChargeField(ptch,Xptch,Y)
    M=length(ptch);
    np=size(Y,1);
    Rptch=zeros(np,M);
    E=zeros(np,3,M);
    for i=1:M
        % r=Xptch(i,:)-Y;
        r=Y-Xptch(i,:);
        Rptch(:,i)=sqrt(sum(r.^2,2));
        E(:,:,i)=ptch(i)/4/pi./(Rptch(:,i).^3).*r;
    end
    E=sum(E,3);
end
