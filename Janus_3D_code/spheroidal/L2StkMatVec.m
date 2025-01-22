function [Stk_x,Stk_y,Stk_z]=L2StkMatVec(pars,sigma_x,sigma_y,sigma_z,Xeval)
%--------------------------------------------------------------------%
% 
%--------------------------------------------------------------------%

if isempty(pars.u0)
    error("No surface parameter u_0 given")
end
if isempty(pars.a)
    error("No scale (a) given")
end
if isempty(pars.matvec_eta)
    error("No matvec_eta provided.");
end

u0=pars.u0;
a=pars.a;
ns=size(pars.centers,1);
if_oblate=pars.oblate;
p=pars.p;
if p==0
    error("p=0.");
end
np=2*p*(p+1);
if length(u0)==1
    u0=u0*ones(1,ns); 
end
if length(a)==1
    a=a*ones(1,ns); 
end

if nargin==4
    % Self evaluation of SL_Stk
    [Stk_x,Stk_y,Stk_z]=L2Stk([],pars,sigma_x,sigma_y,sigma_z,ns);

    % TODO: all to self / self to all.

else

    nt=size(Xeval,1);

    [sep,Xt]=pars.separate_targets(Xeval);

    % -----------
    % fprintf("\n number of separated targets = %d",sum(sep));
    %-------------------------------------------
    
    % Near-field points calculated via spectral method

    X_spectral = cell(1,ns);
    for i=1:ns
        X_spectral{i}=Xt(sep(:,i)==0,:,i);
    end
    
    [Stkspectral_x_cell,Stkspectral_y_cell,Stkspectral_z_cell]=L2Stk(X_spectral,pars,sigma_x,sigma_y,sigma_z,ns);

    % Far evaluation with smooth quadrature.

    Stk_x=zeros(nt,1); Stk_y=Stk_x; Stk_z=Stk_x;

    % Calculate effect from each particle on each target point

    for i=1:ns
        % Add contributions from spectral method
        Stk_x(sep(:,i)==0,:) = Stk_x(sep(:,i)==0,:) + Stkspectral_x_cell{i}; 
        Stk_y(sep(:,i)==0,:) = Stk_y(sep(:,i)==0,:) + Stkspectral_y_cell{i}; 
        Stk_z(sep(:,i)==0,:) = Stk_z(sep(:,i)==0,:) + Stkspectral_z_cell{i}; 
    
        % Get targets for smooth quadrature
        X_smooth = Xt(sep(:,i)==1,:,i);
        nt_smooth=size(X_smooth,1);

        if nt_smooth>0
            if ~if_oblate(i)
                Xself=prolate_spheroid_shape(p,u0(i),a(i));
            else
                Xself=oblate_spheroid_shape(p,u0(i),a(i));
            end
            Sns=SurfaceSph(Xself);

            KEparams = Kernel_Eval_parameters('SL_Stk_3D',0,1,1,1,1e-8,2,400,1);
            KEparams.dim = 3; KEparams.mu=1;
            [~, gwt]=g_grid(p+1);
            wt = pi/p*repmat(gwt', 2*p, 1)./sin(gl_grid(p));
            wt = wt(:);
            Wns = Sns.geoProp.W; Wns= Wns.*wt;
            Wv = repmat(Wns,1,3)'; Wv=Wv(:);
            Xv = reshape(repmat(Xself,1,3)',3,[])';
            Xtrg_ii = reshape(repmat(X_smooth,1,3)',3,[])';
            KEparams.X = Xv;
            KEparams.W2 = Wv.';
            KEparams.cj=repmat((1:3)',np,1);
            KEparams.ci=repmat((1:3)',nt_smooth,1);

            LP_Kernel = Kernel_Eval(Xtrg_ii,Xv,KEparams);

            sig = reshape([sigma_x(:,:,i),sigma_y(:,:,i),sigma_z(:,:,i)].',[],1);
            LP = LP_Kernel*sig;
    
            LP=reshape(LP,3,[]).';
            Stk_x(sep(:,i)==1,:) = Stk_x(sep(:,i)==1,:) + LP(:,1);
            Stk_y(sep(:,i)==1,:) = Stk_y(sep(:,i)==1,:) + LP(:,2);
            Stk_z(sep(:,i)==1,:) = Stk_z(sep(:,i)==1,:) + LP(:,3);

       end
    end

end





