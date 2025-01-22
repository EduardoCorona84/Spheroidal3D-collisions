% A comprehensive test script for convergence of Layer potential operator
% spectra compared against SDY(on surface) and Kernel_Eval(off surface). 
% Reference: DLtest, SPtest by L Crowder, DLtest_oblate by T Li

function spectra_conv_test(bie,obl)
    %{
    bie: layer operator to be checked,
        {'DL','SL','SP'}
    obl: boolean of whether the source is an oblate or not.
        {0,1}
    %}
    if nargin==0
        bie='DL';
        obl=0;
    elseif nargin==1
        obl=0;
    end
    pstart=2; pend=22; 
    % pstart=8;pend=8;
    parr=pstart:2:pend; 
    % Create sources spheroid
    uu=2/sqrt(3);
    % uu=5/3;
    if ~obl
        a=1./uu;
        obl_str='Prolate';
        % Xself=prolate_spheroid_shape(pend,uu,a);
    else
        a = 1./sqrt(1+uu.^2);
        obl_str='Oblate';
        % Xself=oblate_spheroid_shape(pend,uu,a);
    end
    f = @(u,v) exp(sin(u).^4.*cos(u).*cos(4*v));

    % Test 1: one target prolate spheroid
    ns=1;
    uutrg=1.1.*uu;
    atrg=1./uutrg;
    pars=SpheroidalParameters;
    pars.u0 = uu;
    pars.a = a;
    pars.oblate=obl; 
    Xtrg=prolate_spheroid_shape(pend,uutrg,atrg);
    dists=[2,3,4];
    Xtrg_dist=Xtrg + dists.*ones(size(Xtrg));
    if strcmp(bie,'SP')
        % Strg=SurfaceSph(Xtrg_dist); 
        % nu_vec=reshape(Strg.geoProp.nor.to_array,[],3);
        % Nu_cell=nu_vec;

        nu_x=repmat([1,0,0],size(Xtrg_dist,1),1,1);
        nu_y=repmat([0,1,0],size(Xtrg_dist,1),1,1);
        nu_z=repmat([0,0,1],size(Xtrg_dist,1),1,1);
        Nu_cell=nu_x;
    end

    % % Test 2: 5 target prolate spheroids, disregarding overlaps among
    % % targets, each at a different distance and different direction out.
    % ns=5;
    % obl_list=randi(2,ns,1)-1;
    % uulist=linspace(1.1,1.5,ns).*uu;
    % alist=zeros(size(uulist));
    % for trg_ind=1:ns
    %     if obl_list(trg_ind)
    %         alist(trg_ind)=1./sqrt(1+uulist(trg_ind)^2);
    %     else
    %         alist(trg_ind)=1./uulist(trg_ind);
    %     end
    % end
    % % alist=1./uulist;
    % dists=linspace(5,8,ns)';
    % sph_center=[dists,2.*(rand(ns,1)-0.5),2.*pi.*rand(ns,1)];
    % centers=spheroidal2cart(sph_center,a,obl);
    % Xtrgs=cell(1,ns);
    % Nu_cell=cell(1,ns);
    % for tind=1:ns
    %     if obl_list(tind)
    %         Xtrg_base=oblate_spheroid_shape(pend,uulist(tind),alist(tind));
    %     else
    %         Xtrg_base=prolate_spheroid_shape(pend,uulist(tind),alist(tind));
    %     end
    %     Xtrg_dist=Xtrg_base + centers(tind,:);
    %     Xtrgs{tind}=Xtrg_dist;
    %     if strcmp(bie,'SP')
    %         Strg=SurfaceSph(Xtrg_dist); 
    %         Nu_cell{tind} = reshape(Strg.geoProp.nor.to_array,[],3);
    %     end
    % end
    % pars=SpheroidalParameters;
    % pars.u0=uu; pars.a=a;
    % pars.oblate=obl; 
    
    %% on surface
    % TEST1: Compare with SPY
    % %{
    errnorm = zeros(length(parr),1);
    for l=1:length(parr)
        p=parr(l);
        [u,v]=gl_grid(p);
        Y = f(u,v); 
        pars.sigma= Y;
    
        % using kernel matrix 
        [~,dSY] = SDY(p,bie,obl,Y);
        % using spectra
        % test target point on surface eval %%%%%%%%%%%%%%%
        [Xi,~]=pars.get_X();
        % Note: only works for SL, DL, when one spheroid is present.
        %%%%%%%%%%%%%%%%%%%%%%%%%

        switch bie
            case 'DL'
                myDY=spheroidalDL(pars,Xi);
            case 'SL'
                myDY=spheroidalSL(pars,Xi);
                % myDY2=spheroidalSL(pars);
            case 'SP'
                myDY=spheroidalSP(pars);
        end
        errnorm(l) = norm(abs(myDY - dSY));
    end
    plot_name=append('On Surface Convergence of ',obl_str,' Source ',bie,' Spectra');
    figure;
    semilogy(parr,errnorm);
    title(plot_name);
    % %}

    % TEST2: Compare with higest p
    %{
    errnorm = zeros(length(parr),1);
    errnorm_y = zeros(length(parr),1);
    errnorm_z = zeros(length(parr),1);
    for l=length(parr):-1:1
        p=parr(l);
        [u,v]=gl_grid(p);
        Y = f(u,v); 
        pars.sigma= Y; pars.get_shc();
        if p==pend
            [Xi,~]=pars.get_X();
        end

        switch bie
            case 'DL'
                myDY=spheroidalDL(pars,Xi);
            case 'SL'
                myDY=spheroidalSL(pars,Xi);
                % myDY2=spheroidalSL(pars);
            case 'SP'
                nu_x=repmat([1,0,0],size(u,1),1,ns);
                nu_y=repmat([0,1,0],size(u,1),1,ns);
                nu_z=repmat([0,0,1],size(u,1),1,ns);
                [myDY,myDYy,myDYz]=spheroidalSP(pars,[],nu_x,nu_y,nu_z);
        end
        if p==pend
            dSY=myDY;
            if strcmp(bie,'SP')
                dSYy=myDYy;
                dSYz=myDYz;
            end
        end
        % TODO: dimensions don't match. have to match grid points.
        errnorm(l) = norm(abs(myDY - dSY));
        if strcmp(bie,'SP')
            errnorm_y(l) = norm(abs(myDYy-dSYy));
            errnorm_z(l) = norm(abs(myDYz-dSYz));
        end
    end

    figure;
    semilogy(parr(1:end-1),errnorm(1:end-1));
    if strcmp(bie,'SP')
        hold on;
        semilogy(parr(1:end-1),errnorm_y(1:end-1));
        semilogy(parr(1:end-1),errnorm_z(1:end-1));
        hold off;
        legend('dim 1','dim 2','dim 3');
    end
    %}

    %% off surface
    %{
    [up,vp]=gl_grid(pend);
    Yp=f(up,vp);
    
    % TEST1: Compare with Kernel_Eval
    %{
    switch bie
        case 'DL'
            pot='DL_L_3D';
        case 'SL'
            pot='SL_L_3D';
        case 'SP'
            pot='dSL_L_3D';
    end
    KEparams = Kernel_Eval_parameters(pot,0,1,1,1,1e-8,2,400,1);
    KEparams.dim = 3;
    KEparams.X = Xself;
    Sns = SurfaceSph(Xself);
    [~, gwt]=g_grid(pend+1);
    wt = pi/pend*repmat(gwt', 2*pend, 1)./sin(gl_grid(pend));
    wt = wt(:);
    Wns = Sns.geoProp.W; Wns= Wns.*wt;
    KEparams.W2 = Wns.';

    % 1 spheroid
    err_norm=zeros(length(parr),1); err_norm_y=err_norm; err_norm_z=err_norm;
    for k=1:length(parr)
        p=parr(k);
        [u,v]=gl_grid(p);
        Y = f(u,v); 
        pars.sigma=Y;pars.get_shc();

        switch bie
            case 'DL'
                myDY=spheroidalDL(pars,Xtrg_dist);
                Nrns = reshape(Sns.geoProp.nor.to_array,[],3);
                KEparams.nor = Nrns; 
            case 'SL'
                myDY=spheroidalSL(pars,Xtrg_dist);
            case 'SP'
                % myDY=spheroidalSP(pars,Xtrg_dist,Nu_cell);
                KEparams.nor=Nu_cell;
                [myDY,SP2,SP3]=spheroidalSP(pars,Xtrg_dist,Nu_cell,nu_y,nu_z);
        end
        Smat = Kernel_Eval(Xtrg_dist,Xself,KEparams); 
        mat_mulY=Smat * Yp;
        err_norm(k)=norm(abs(myDY-mat_mulY));

        if strcmp(bie,'SP')
            KEparams.nor=nu_y;
            Smat = Kernel_Eval(Xtrg_dist,Xself,KEparams); 
            mat_mulY=Smat * Yp;
            err_norm_y(k)=norm(abs(SP2-mat_mulY));
    
            KEparams.nor=nu_z;
            Smat = Kernel_Eval(Xtrg_dist,Xself,KEparams); 
            mat_mulY=Smat * Yp;
            err_norm_z(k)=norm(abs(SP3-mat_mulY));
        end
    end
    

    % % 5 spheroids
    % MatVeccell=cell(1,ns);
    % err_norm=zeros(length(parr),ns);
    % for tind=1:ns
    %     Xtrg_t=Xtrgs{tind};
    %     switch bie
    %         case 'DL'
    %             Nrns = reshape(Sns.geoProp.nor.to_array,[],3);
    %             KEparams.nor = Nrns; 
    %         case 'SP'
    %             nu_cart=Nu_cell{tind};
    %             KEparams.nor=nu_cart;
    %     end
    %     Smat = Kernel_Eval(Xtrg_t,Xself,KEparams); 
    %     mat_mulY=Smat * Yp;
    %     MatVeccell{tind}=mat_mulY;
    % 
    %     for k=1:length(parr)
    %         p=parr(k);
    %         [u,v]=gl_grid(p);
    %         Y = f(u,v); 
    %         pars.sigma=Y;
    %         switch bie
    %             case 'DL'
    %                 myDY=spheroidalDL(pars,Xtrg_t);
    %             case 'SL'
    %                 myDY=spheroidalSL(pars,Xtrg_t);
    %             case 'SP'
    %                 myDY=spheroidalSP(pars,Xtrg_t,nu_cart);
    %         end
    %         err_norm(k,tind)=norm(abs(myDY-mat_mulY));
    %     end
    % end

    % plot_name=append('Off Surface Convergence of ',obl_str,' Source ',bie,' Spectra');
    % figure;
    % for tind=1:ns
    %     err_p=err_norm(:,tind);
    %     semilogy(parr,err_p); hold on;
    % end
    % hold off;
    % xlabel("p"); ylabel("norm error");
    % title(plot_name)
    % 
    % if strcmp(bie,'SP')
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     figure;
    %     for tind=1:ns
    %         err_p=err_norm_y(:,tind);
    %         semilogy(parr,err_p); hold on;
    %     end
    %     hold off;
    % 
    %     figure;
    %     for tind=1:ns
    %         err_p=err_norm_z(:,tind);
    %         semilogy(parr,err_p); hold on;
    %     end
    %     hold off;
    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % end

    %}


    %}
end
