function [Stk_XX,Stk_YY,Stk_ZZ]=test_stokes(ns,p,plot)
    % Curvature Stokian flow over one/large system of prolates/mixes

    if nargin<3
        plot=0;
    end

    % NOTE: this ./data is the data folder inside spheroidal, not on same
    % level of it.
    setupName='./data/stokes_'+string(ns)+'_ptcls_p'+string(p)+'_setup.txt';
    stksxName='./data/stokes_'+string(ns)+'_ptcls_p'+string(p)+'x.txt';
    stksyName='./data/stokes_'+string(ns)+'_ptcls_p'+string(p)+'y.txt';
    stkszName='./data/stokes_'+string(ns)+'_ptcls_p'+string(p)+'z.txt';

    np=2*p*(p+1);

    rFlag=exist(setupName,'file');
    if rFlag
        % ----- read data from file ---------------------------------------
        [p,ns,centers,u0s,as,oblates,Rmats] = readSystem(setupName);
        pars=SpheroidalParameters;
        pars.matvec_eta = 10;
        pars.isReal=0;
        pars.u0=u0s; 
        pars.a = as;
        pars.oblate=oblates;
        pars.centers=centers; 
        pars.Rmat=Rmats;
        % dbox=pars.a(end)*4/sqrt(7);
        dbox=pars.a(end)*5/4;

    else
        % ----- Set up system of spheroids --------------------------------
        if ns==3
            ns=ns+1;
            % u0=[1.1 1.2 1.3 3/sqrt(7)];
            u0=[1.1 1.2 1.3 5/4];
            % ns=ns;
            % u0=[1.1 1.2 1.3];
            pars=SpheroidalParameters;
            % pars.matvec_eta = 10;
            pars.matvec_eta = 2;
            pars.isReal=0;
            pars.u0=u0;
            pars.p=p;
            alist = 1./u0;
            % alist(end) = 8./sqrt(u0(end)^2+1);
            alist(end) = 8./u0(end);
            pars.a=alist;
            % pars.oblate = [0 0 0 1];
            pars.oblate=[0 0 0 0];
        
            pars.centers = [0 0 0; 4 0 0; 1.6 1.6 1.6; 0 0 0];
            thetas = [0 pi/10 5*pi/3 pi/2];
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
            dbox=8;
        
        elseif ns==1
            ns=ns;
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
            dbox=1.5;
        
            pars.Rmat=eye(3,3);
        
        else
            [u0,a,oblate,centers,Rmats,ns_actual,dbox]=setup_large(ns,p);
            ns=ns_actual+1;
            pars=SpheroidalParameters;
            pars.matvec_eta = 10;
            pars.isReal=0;
            pars.u0=[u0 5/4];
            pars.a = [a dbox*4/5];
            pars.oblate = [oblate 0];
            pars.Rmat=cat(3,Rmats,[0 0 1;0 1 0;-1 0 0]);

            % pars.u0=[u0 3/sqrt(7)]; 
            % pars.a = [a dbox*sqrt(7)/4];
            % pars.oblate=[oblate 1];
            % pars.Rmat=cat(3,Rmats,eye(3));

            pars.centers=[centers; 0 0 0]; 
        end
        %------------------------------------------------------------------

        % % ------- Writing setup to file -----------------------------------
        % fID=fopen(setupName,'w');
        % fprintf(fID,'p=%d\n',p); % changed: store p value as well.
        % fprintf(fID,"ns total=%d\n",ns); % changed: now stores total number of spheroids (including outer shell)
        % 
        % fprintf(fID,"centers:x,y,z\n");
        % fprintf(fID,'%f, %f, %f\r\n',pars.centers');
        % fprintf(fID,"\n ");
        % 
        % fprintf(fID,'u0s:\n');
        % fprintf(fID,'%f\r\n',pars.u0);
        % fprintf(fID,"\n ");
        % 
        % fprintf(fID,"a's:\n");
        % fprintf(fID,'%f\r\n',pars.a);
        % fprintf(fID,"\n ");
        % 
        % fprintf(fID,'oblate:\n');
        % fprintf(fID,'%i\r\n',pars.oblate);
        % fprintf(fID,"\n ");
        % 
        % fprintf(fID,'Rmat: x,y,z\n');
        % for i=1:ns
        %     Mat_temp = pars.Rmat(:,:,i);
        %     fprintf(fID,'%f, %f, %f\r\n',Mat_temp');
        %     fprintf(fID,'\n');
        % end
        % fclose(fID);


    end
    % ---------------------------------------------------------------------

    % ------------- For curvature flow, sigma=H*n -------------------------
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
        sigma_x(:,:,i)=K.*Nri(:,1); % Nri outward normal on spheroids.
        sigma_y(:,:,i)=K.*Nri(:,2);
        sigma_z(:,:,i)=K.*Nri(:,3);
        sigma_norm=sqrt((K.*Nri(:,1)).^2+(K.*Nri(:,2)).^2+(K.*Nri(:,3)).^2);
        % if plot
        %     if i~=ns
        %         plotb([Yi(:,1);Yi(:,2);Yi(:,3)],sigma_norm); hold on;
        %     else
        %         Yout=prolate_spheroid_shape(24,pars.u0(end),pars.a(end));
        %         Yout=(pars.Rmat(:,:,end)*Yout')';
        % 
        %         % Not the same as 'ellipseX'! 
        %         % scatter3(Yout(:,1),Yout(:,2),Yout(:,3),'bo'); hold on;
        %         % YS=SurfaceSph(shape_gallery(p*2,'ellipseX'));
        %         % Yout1=8.*reshape(YS.cart.to_array(),[],3);
        %         % scatter3(Yout1(:,1),Yout1(:,2),Yout1(:,3),'r*');
        % 
        %         sigma_norm=interpsh(sigma_norm,24);
        %         outer=plotb([Yout(:,1);Yout(:,2);Yout(:,3)],sigma_norm); hold on;
        %         set(outer,'facealpha',0.2);
        %         axis equal
        %         axis off;
        % 
        %         cb = colorbar; 
        %         set(cb,'position',[.15 .2 .05 .3])
        %         % hold off;
        % 
        %     end
        % 
        % end

    end


    
    % Second plot: cross section at z=0
    if plot
        figure;
        for i=1:ns-1
            Yi=Y((i-1)*np+1:i*np,:);
            zi=Yi(:,3);
            % pos=zi(1)>1e-12; % whether first zi is above or below z=0
            pos=(zi(1)-2)>1e-12;
            zind=2;

            while zind<length(zi) % go through zi
                if pos~=((zi(zind)-2>1e-12))
                    if pars.oblate(i)
                        Ysmooth=oblate_spheroid_shape(24,pars.u0(i),pars.a(i));
                    else
                        Ysmooth=prolate_spheroid_shape(24,pars.u0(i),pars.a(i));
                    end
                    ci=pars.centers(i,:)';
                    Ysmooth =( pars.Rmat(:,:,i) * Ysmooth' + repmat(ci,1,2*24*25) )';
                    plotb([Ysmooth(:,1);Ysmooth(:,2);Ysmooth(:,3)]); hold on; 
                    break
                end
                zind = zind + 1;
            end

        end   
    end

    %{
    % ----- Set one evaluation point to check close vs far evaluation -----
    Xeval=[-6.2,-7.12,2.95];
    matvec_to_test=[6,8,10,12];
    pars.matvec_eta=14;
    [Stk_x_eval,Stk_y_eval,Stk_z_eval]=L2Stk(Xeval,pars,sigma_x,sigma_y,sigma_z,ns,theta,phi);
    err_x=zeros(1,4); err_y=err_x; err_z=err_x;
    for ii=1:4
        pars.matvec_eta=matvec_to_test(ii);
        [Stk_x,Stk_y,Stk_z]=L2Stk(Xeval,pars,sigma_x,sigma_y,sigma_z,ns,theta,phi);
        err_x(ii)=norm(norm(Stk_x_eval-Stk_x));
        err_y(ii)=norm(norm(Stk_y_eval-Stk_y));
        err_z(ii)=norm(norm(Stk_z_eval-Stk_z));
    end
    display(err_x); display(err_y); display(err_z);
    % ---------------------------------------------------------------------
    %}


    % %{
    % ----- Solve for Stokes ----------------------------------------------
    % ----- Set up domain points for velocity field -----------------------

    fprintf("\n size of outer shell, dbox=%f\n",dbox);
    DX=(-dbox+0.02:0.106:dbox+0.02)'; 
    DY=(-dbox*3/5+0.02:0.106:dbox*3/5+0.02)'; DZ=[2];
    [XX,YY,ZZ]=meshgrid(DX,DY,DZ);
    Xbox=[XX(:),YY(:),ZZ(:)];

    
    is_clear_bound=check_clear(pars.centers(end,:),pars.u0(end),pars.a(end),pars.oblate(end),pars.Rmat(:,:,end),Xbox);
    Xgrid=Xbox(is_clear_bound==0,:);
    fprintf("\n size of grid inside shell, should be all of Xbox: %d\n",size(Xgrid,1));

    % only look at points that are also outside the interior spheroids.
    clear_mat=ones(size(Xgrid,1),1);
    for i=1:ns-1 % check overlap with every particle inside outer constraint.
        is_clear=check_clear(pars.centers(i,:),pars.u0(i),pars.a(i),pars.oblate(i),pars.Rmat(:,:,i),Xgrid);
        clear_mat=clear_mat.*is_clear; % if 0 (not clear) for one spheroid i, will be 0 at the end.
    end
    Xeval=Xgrid(clear_mat==1,:);

    fprintf("\n size of eval points: %d\n", size(Xeval,1));

    % ---------------------------------------------------------------------

    rFlag=exist(stksxName,'file');
    % rFlag=0;
    if rFlag
        fprintf("\n opening existing file\n");
        % read Stokes matrices from file as well
        fID=fopen(stksxName);
        % Stk_XX=fscanf(fID,'%f %f %f',[3,Inf])';
        Stk_XX_1d=fscanf(fID,'%f',[1,Inf])';
        fclose(fID);

        fID=fopen(stksyName);
        % Stk_YY=fscanf(fID,'%f %f %f',[3,Inf])';
        Stk_YY_1d=fscanf(fID,'%f',[1,Inf])';
        fclose(fID);

        fID=fopen(stkszName);
        % Stk_ZZ=fscanf(fID,'%f %f %f',[3,Inf])';
        Stk_ZZ_1d=fscanf(fID,'%f',[1,Inf])';
        fclose(fID);

        Stk_XX=reshape(Stk_XX_1d,size(XX)); 
        Stk_YY=reshape(Stk_YY_1d,size(YY));
        Stk_ZZ=reshape(Stk_ZZ_1d,size(ZZ));

    else
        % --------- Stokes Layer potential Evaluation ---------------------
        Stk_x=zeros(size(Xgrid,1),1); Stk_y=Stk_x; Stk_z=Stk_x;
        tic;
        % [Stk_x_eval,Stk_y_eval,Stk_z_eval]=L2Stk(Xeval,pars,sigma_x,sigma_y,sigma_z,ns);
        [Stk_x_eval,Stk_y_eval,Stk_z_eval]=L2StkMatVec(pars,sigma_x,sigma_y,sigma_z,Xeval);
        
        tend=toc;
        fprintf("\n time to perform Stokes calculation for p=%d: %d seconds.\n", p, tend);

        Stk_x(clear_mat==1,:)=real(sum(Stk_x_eval,3));
        Stk_y(clear_mat==1,:)=real(sum(Stk_y_eval,3));
        Stk_z(clear_mat==1,:)=real(sum(Stk_z_eval,3));

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
        % -----------------------------------------------------------------
        
        % ------ Write Stokes matrices to file. ---------------------------
        % reshape to 1 x nd^3 array to store
        Stk_XX_1d=reshape(Stk_XX,[],1);
        Stk_YY_1d=reshape(Stk_YY,[],1);
        Stk_ZZ_1d=reshape(Stk_ZZ,[],1);

        fID=fopen(stksxName,'w');
        % fprintf(fID,'%f %f %f\r\n',Stk_XX);
        fprintf(fID,'%f\r\n',Stk_XX_1d);
        fclose(fID);

        fID=fopen(stksyName,'w');
        % fprintf(fID,'%f %f %f\r\n',Stk_YY);
        fprintf(fID,'%f\r\n',Stk_YY_1d);
        fclose(fID);

        fID=fopen(stkszName,'w');
        % fprintf(fID,'%f %f %f\r\n',Stk_ZZ);
        fprintf(fID,'%f\r\n',Stk_ZZ_1d);
        fclose(fID);

    end
    
    if plot

        % x0=linspace(min(Xbox(:,1)),max(Xbox(:,1)),100); x1=linspace(min(Xbox(:,2)),max(Xbox(:,2)),100);x2=linspace(min(Xbox(:,3)),max(Xbox(:,3)),100); 
        % [XX0,XX1,XX2]=meshgrid(x0,x1,x2);
        % UU0=griddata(XX,YY,ZZ,Stk_XX,XX0,XX1,XX2);
        % UU1=griddata(XX,YY,ZZ,Stk_YY,XX0,XX1,XX2);
        % UU2=griddata(XX,YY,ZZ,Stk_ZZ,XX0,XX1,XX2);
        % UU0(isnan(UU0))=0; UU1(isnan(UU1))=0; UU2(isnan(UU2))=0;
        % streamslice(XX0,XX1,XX2,UU0,UU1,UU2,[],[],XX2(1,1,50),10);

        % slice=33;
        % ZZ_slice=ZZ(:,:,slice); zval=ZZ_slice(1,1);
        % VV_slice=VV(slice,:,:); vval=VV_slice(1,1);
        % title("Flow line for z="+string(zval));
        zval=2;
        % HS = streamslice(XX,YY,ZZ,Stk_XX,Stk_YY,Stk_ZZ,[],[],zval,10);
        HS = streamslice(XX,YY,Stk_XX,Stk_YY,10); hold on;
        set(HS, 'LineWidth', 0.1, 'Color', [128 128 128]/255); %axis off;  hold on
        axis equal;
        axis off;
        % outer shell outline
        tt=linspace(0,2*pi+2*pi/100,100);
        rout=sqrt(dbox^2-25*zval^2/9);
        xouter=rout.*cos(tt);
        youter=3/5*rout.*sin(tt);
        zouter=zval.*ones(size(xouter));
        plot3(xouter,youter,zouter,'k'); 
        % scatter3(Xeval(:,1),Xeval(:,2),Xeval(:,3));
        % quiver3(XX(:,:,33),YY(:,:,33),ZZ(:,:,33),Stk_XX(:,:,33),Stk_YY(:,:,33),Stk_ZZ(:,:,33))
        view(2);
        hold off;
        
        % steps=floor(length(DZ)/12);
        % if steps~=0
            % --- 2D streamline plotting ----------------------------------
            % figure;
            % slices=floor((length(DZ)-12*steps)/2)+steps.*(1:12);
            % slices=[10,11,12,13,25,26,27,28,34,35,36,37];
            % nslice=length(slices);
            % for s=1:nslice
                % subplot(3,4,s)
                % ZZ_slice=ZZ(:,:,slices(s));
                % zval=ZZ_slice(1,1);
                % title("Flow line for z="+string(zval));
                % % HS = streamline(stream2(XX(:,:,slices(s)), YY(:,:,slices(s)), Stk_XX(:,:,slices(s)), Stk_YY(:,:,slices(s)), XX(:,:,slices(s)), YY(:,:,slices(s))));
                % % HS = streamline(stream2(XX(:,:,slice), YY(:,:,slice), Stk_XX(:,:,slice), Stk_YY(:,:,slice), XX(:,:,slice), YY(:,:,slice)));
                % HS = streamslice(XX,YY,ZZ,Stk_XX,Stk_YY,Stk_ZZ,[],[],zval,2.5);
                % set(HS, 'LineWidth', 0.1, 'Color', 'r'); %axis off;  hold on
                % % set(HS, 'LineWidth', 0.02); %axis off;  hold on
                % axis equal;
                % % set(HS, 'LineWidth',1); 
                % xlabel("x");ylabel("y");
                % % plot([XX(:,1,slice); XX(1,1,slice)], [YY(:,1,slice); YY(:,1,slice)],'r');
                % hold off;
            % end
            % hold off;

            % pic_name='./data/stokes_'+string(ns)+'_ptcls_streamline.png';
            % fstream=gcf;
            % exportgraphics(fstream,pic_name,'Resolution',300);

            % streamslice()

            % % --- 3D streamline plotting -----------------------------
            % figure;
            % HS=streamline(XX,YY,ZZ,Stk_XX,Stk_YY,Stk_ZZ,XX,YY,zeros(size(XX)));
            % % HS=streamline(stream3(XX,YY,ZZ,Stk_XX,Stk_YY,Stk_ZZ,XX,YY,ZZ));
            % set(HS,'LineWidth',1);

        % end


    end
    % %}


    function is_clear = check_clear(cntr,u0,a,ob,rot_mat,xtrg)
        % check if points xtrg are inside a given spheroid
    
        % input:
        %       cntr: 1 x 3 position of spheroid.
        %       u0, a, ob: parameters of the spheroid i.
        %       xtrg: n x 3 points to determine.
    
        % output: 
        %       is_clear: n x 1 boolean of whether each point is "clear", i.e.
        %       outside the given spheroid.

        xtrg_rot=(rot_mat'*(xtrg-cntr)')'+cntr;
    
        xdis=xtrg_rot(:,1)-cntr(1);
        ydis=xtrg_rot(:,2)-cntr(2);
        zdis=xtrg_rot(:,3)-cntr(3);
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

end



