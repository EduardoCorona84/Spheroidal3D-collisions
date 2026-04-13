function A = Kernel_Eval(X1,X2,params)
%
% This file is part of HSSDirectSolver
% Copyright (C) 2013 Eduardo Corona, Denis Zorin, Per Gunnar Martinsson
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%     FUNCTION CALL:
%         A = Kernel_Eval(X1,X2,params)
%
%     top - for nonsym multiplication by b and c
%     DESCRIPTION:
%         This function takes 2D point arrays X1 and X2, and returns the matrix
%         A = K[X1,X2].
%
%     INPUT:
%         params (see HSS_tree_parameters.m)
%         params.flagpot - {'SL_H_2D', 'SL_H_3D', 'SL_L_2D', 'SL_L_3D',
%         'SL_Y_2D', 'SL_Y_3D'}
%

% Modify params struct accordingly
[params,X1,X2] = set_params(params,X1,X2); 

if ~isfield(params,'W2')
    h = params.h; 
    wh = params.h^2;  
else
    wh = repmat(params.W2,size(X1,1),1);   
end

%display('size of weight is:');
%display(size(wh));

%old code
%{
%X grid
[Y_g1,  X_g1  ] = meshgrid(X2(:,1), X1(:,1));
%}


if params.dim == 1
   den = abs(X_g1 - Y_g1); 
   switch params.flag_pot
       case 'SL_L_2D'
          wh = (1/2/pi)*wh; 
          A = (den==0) + 0.5*wh.*log((den + (den==0))); 
   
   end

% 2D (plane or curves on plane)
elseif params.dim == 2
    % Y grid 
    [Y_g2,  X_g2  ] = meshgrid(X2(:,2), X1(:,2)); 
     % Determine functions b(x) and c(y)
    [b,c] = LOCAL_get_bc({X_g1,X_g2},{Y_g1,Y_g2},params,4);
    % den = ||X-Y||^2
     den = (X_g1 - Y_g1).^2 + (X_g2 - Y_g2).^2;
    
    switch params.flag_pot
     
    % 2D Helmholtz
    case 'SL_H_2D'
        kh = params.kh; 
        C = besselh(0,kh);
        if params.order == 4
            w = 1 + kh^2*(0.25*1i)*params.dr_weights; 
            A = w*(den==0) + wh*kh^2*(0.25*1i)*b.*(besselh(0,kh*sqrt(den + (den==0))) - C*(den==0)).*c;
        else
            % Have to add q-diagonal correction here. 
            A = (den==0) + wh*kh^2*(0.25*1i)*b.*(besselh(0,kh*sqrt(den + (den==0))) - C*(den==0)).*c;
        end
    % 3D Helmholtz    
    case 'SL_H_3D'
        kh = params.kh; 
        C = exp(1i*kh);
        A = (den==0) + wh*kh*(exp(1i*kh*sqrt(den + (den==0)))./sqrt(den + (den==0)) - C*(den==0)); 
    % 2D Laplace    
    case 'SL_L_2D'
        wh = (1/2/pi)*wh; 
        A = (den==0) + 0.5*wh.*b.*log((den + (den==0))).*c;
    % 3D Laplace    
    case 'SL_L_3D'
        wh = (1/4/pi)*wh; 
        A = (den==0) + wh.*(1./sqrt(den + (den==0)) - (den==0)); 
    % Differential Operator Laplace (5pt stencil)
  
    case 'DO_L_2D_5' 
        den = (1/h)*sqrt(den);     
        A = -4*(den==0)+1*(den==1); 
    % Differential Operator Laplace (9pt stencil)    
    case 'DO_L_2D_9'
        den = (1/h)*sqrt(den);      
        A = -(20/6)*(den==0)+(4/6)*(den==1)+(1/6)*(den==sqrt(2));       
    % 2D Yukawa    
    case 'SL_Y_2D'
        C = besselk(0,params.kh);  
        A = (den==0) + wh*(0.5/pi)*b.*(besselk(0,params.kh*(den + (den==0))) - C*(den==0)).*c; 
    %3D Yukawa    
    case 'SL_Y_3D'
        A = (den==0) + wh*b.*(exp(-params.kh*den)./sqrt(den + (den==0)) - (den==0)).*c;  
    case 'Distance'
        A = sqrt(den); 
    case 'zeros'
         A = zeros(size(den)); 
    case 'Fun_2D'
        fun =params.fun; 
        A  =wh.*fun(den);  
    case 'Kernel_Name'
        % User-defined Kernel template K(x,y) = I + b(X)K(|X-Y|)c(Y) 
        A = (den==0) + 0.5*wh.*b(X_g1,X_g2).*K(den).*c(Y_g1,Y_g2);
    end
 % 3D (surfaces and volume)   
else

    %old code
    
    % Y grid 
    %{
    [Y_g1,  X_g1  ] = meshgrid(X2(:,1), X1(:,1));
    [Y_g2,  X_g2  ] = meshgrid(X2(:,2), X1(:,2)); 
    % Z grid
    [Y_g3,X_g3] = meshgrid(X2(:,3),X1(:,3));
    
    % den = ||X-Y||^2
    d1 = (X_g1 - Y_g1); d2 = (X_g2 - Y_g2); d3 = (X_g3 - Y_g3); 
    den = d1.^2 + d2.^2 + d3.^2;
    %}


     den = (X1(:,1) - X2(:,1).').^2 + (X1(:,2) - X2(:,2).').^2 + (X1(:,3) - X2(:,3).').^2;

     switch params.flag_pot
     case 'SL_L_3D'  
        % Determine functions b(x) and c(y) 
        wh = (1/(4*pi))*wh; 
        a = params.a; 
        % Single layer Laplace 
   %     display('size of density is:')
  %      display(size(den));
        A = a*(den==0) + wh.*(1./sqrt(den + (den==0)) - (den==0)); 
     case 'DL_L_3D'  
        % Double layer Laplace
        a = params.a;
        wh = (1/(4*pi))*wh;
        N1 = repmat(params.nor(:,1).',size(X1,1),1);
        N2 = repmat(params.nor(:,2).',size(X1,1),1);
        N3 = repmat(params.nor(:,3).',size(X1,1),1);
        NdotR = d1.*N1+d2.*N2+d3.*N3;    
        A = a*(den==0) + wh.*NdotR.*(1./(den+(den==0)).^(3/2) - (den==0));   
     case 'dSL_L_3D'
        % dS/dn of Laplace
        wh = -(1/(4*pi))*wh;      
        a = params.a;
        N1 = repmat(params.nor(:,1),1,size(X2,1));
        N2 = repmat(params.nor(:,2),1,size(X2,1));
        N3 = repmat(params.nor(:,3),1,size(X2,1));
        NdotR = d1.*N1+d2.*N2+d3.*N3;    
        A = a*(den==0) + wh.*NdotR.*(1./(den+(den==0)).^(3/2) - (den==0));    
     case 'dDL_L_3D'
           lambda=0;
          a=params.a;
          wh=-(1/(4*pi))*wh;
          %computes the dot products <r,n_source>,<r,n_target>
        SourceN1 = repmat(params.nor(:,1).',size(X1,1),1);
        SourceN2 = repmat(params.nor(:,2).',size(X1,1),1);
        SourceN3 = repmat(params.nor(:,3).',size(X1,1),1);
        NdotRSource = d1.*SourceN1+d2.*SourceN2+d3.*SourceN3;
        TargN1 = repmat(params.targnor(:,1),1,size(X2,1));
        TargN2 = repmat(params.targnor(:,2),1,size(X2,1));
        TargN3 = repmat(params.targnor(:,3),1,size(X2,1));
        NdotRTarg = d1.*TargN1+d2.*TargN2+d3.*TargN3;  
        dotnorm=TargN1.*SourceN1+TargN2.*SourceN2+SourceN3.*TargN3;
        rbar=sqrt(den); %%%
        expy=exp(-lambda*rbar);  %%%
        A=(lambda^2./rbar + 2*lambda./rbar.^2 + 2./rbar.^3);  %%%
        B=(lambda./rbar + 1./rbar.^2);  %%%
        delsquared=-1./rbar.^3.*NdotRTarg.*NdotRSource + 1./rbar.*dotnorm; %%%
        A1 =(expy./rbar.^2).*(A.*NdotRTarg.*NdotRSource);
        A2=expy.*B.*delsquared;
        A= a*(den==0) + wh.*(A1 - A2) -(den==0);
    case 'SL_H_3D'
            % Determine functions b(x) and c(y)
            [b,c] = LOCAL_get_bc({X_g1,X_g2,X_g3},{Y_g1,Y_g2,Y_g3},params,4); 
            
        % Single layer 3D Helmholtz   
        kh = params.kh; 
        C = exp(1i*kh);
        A = (den==0) + kh*wh.*b.*(exp(1i*kh*sqrt(den + (den==0)))./sqrt(den + (den==0)) - C*(den==0)).*c;  
     % Differential Operator Laplace 3D (7pt stencil)   
     case 'DO_L_3D_7'   
        den = (1/h)*sqrt(den);     
        A = -6*(den==0)+1*(den==1);  
     % Single Layer Stokes 3D   
     case 'SL_Stk_3D' 
    %need to think about a way to not explicity compute the d's
     d1 = X1(:,1) - X2(:,1).';
     d2 = X1(:,2) - X2(:,2).';
     d3 = X1(:,3) - X2(:,3).';
     % Single layer Stokes 
     wh = (1/(8*pi*params.mu))*wh; 
     a = params.a; 
     % Row and column tensor coords
     ci = params.ci; cj = params.cj; 
     [CJ,CI] = meshgrid(cj,ci);     
     I1 = CI==1; I2 = CI==2; I3 = CI==3;     
     J1 = CJ==1; J2 = CJ==2; J3 = CJ==3;       
     
     % 1/r diagonal part 
     A1 = (CI==CJ).*(a*(den==0) + wh.*(1./sqrt(den + (den==0)) - (den==0)));
     
     % r_i*r_j/r^3
     A2 = wh./(sqrt(den+(den==0)).^3);
     A2(I1) = A2(I1).*d1(I1); A2(J1) = A2(J1).*d1(J1);
     A2(I2) = A2(I2).*d2(I2); A2(J2) = A2(J2).*d2(J2);
     A2(I3) = A2(I3).*d3(I3); A2(J3) = A2(J3).*d3(J3);
     
     A = A1+A2; 
     % Double Layer Stokes 3D    
     case 'DL_Stk_3D'  
     % Double layer Stokes 
     wh = (3/(4*pi))*wh; 
     a = params.a; 
     % Row and column tensor coords
     ci = params.ci; cj = params.cj; 
     [CJ,CI] = meshgrid(cj,ci);     
     I1 = CI==1; I2 = CI==2; I3 = CI==3;     
     J1 = CJ==1; J2 = CJ==2; J3 = CJ==3; 
     
     % N(y) dot R
     N1 = repmat(params.nor(:,1).',size(X1,1),1);
     N2 = repmat(params.nor(:,2).',size(X1,1),1);
     N3 = repmat(params.nor(:,3).',size(X1,1),1);
     
     NdotR = wh.*(d1.*N1+d2.*N2+d3.*N3);    
     
     % diagonal part 
     A1 = (CI==CJ).*(a*(den==0));
     
     % r_i*r_j/r^5
     A2 = NdotR./((den+(den==0)).^(5/2));
     A2(I1) = A2(I1).*d1(I1); A2(J1) = A2(J1).*d1(J1);
     A2(I2) = A2(I2).*d2(I2); A2(J2) = A2(J2).*d2(J2);
     A2(I3) = A2(I3).*d3(I3); A2(J3) = A2(J3).*d3(J3);
     
     A = A1+A2; 
     case 'SDL_Stk_3D' 
     % Single layer Stokes 
     whD = (3/(4*pi))*wh; 
     wh = (1/(8*pi*params.mu))*wh; 
     a = params.a; 
     % Row and column tensor coords
     ci = params.ci; cj = params.cj; 
     [CJ,CI] = meshgrid(cj,ci);     
     I1 = CI==1; I2 = CI==2; I3 = CI==3;     
     J1 = CJ==1; J2 = CJ==2; J3 = CJ==3; 
     
     % N(y) dot R
     N1 = repmat(params.nor(:,1).',size(X1,1),1);
     N2 = repmat(params.nor(:,2).',size(X1,1),1);
     N3 = repmat(params.nor(:,3).',size(X1,1),1);
     NdotR = whD.*(d1.*N1+d2.*N2+d3.*N3);
     
     % 1/r diagonal part 
     A1 = (CI==CJ).*(a*(den==0) + wh.*(1./sqrt(den + (den==0)) - (den==0)));
     
     % r_i*r_j/r^3
     A2 = wh./(sqrt(den+(den==0)).^3);
     A2(I1) = A2(I1).*d1(I1); A2(J1) = A2(J1).*d1(J1);
     A2(I2) = A2(I2).*d2(I2); A2(J2) = A2(J2).*d2(J2);
     A2(I3) = A2(I3).*d3(I3); A2(J3) = A2(J3).*d3(J3);
     
     % r_i*r_j/r^5
     A3 = NdotR./((den+(den==0)).^(5/2));
     A3(I1) = A3(I1).*d1(I1); A3(J1) = A3(J1).*d1(J1);
     A3(I2) = A3(I2).*d2(I2); A3(J2) = A3(J2).*d2(J2);
     A3(I3) = A3(I3).*d3(I3); A3(J3) = A3(J3).*d3(J3);
     
     A = A1+A2+A3;     
     case 'dSL_Stk_3D'    
     % dS/dn_x Stokes 
     wh = (-1/(8*pi))*wh; 
     a = params.a; 
     % Row and column tensor coords
     ci = params.ci; cj = params.cj; 
     [CJ,CI] = meshgrid(cj,ci);     
     I1 = CI==1; I2 = CI==2; I3 = CI==3;     
     J1 = CJ==1; J2 = CJ==2; J3 = CJ==3; 
     
     % N(x) dot R
     N1 = wh.*repmat(params.nor(:,1),1,size(X2,1))./(sqrt(den+(den==0)).^3);
     N2 = wh.*repmat(params.nor(:,2),1,size(X2,1))./(sqrt(den+(den==0)).^3);
     N3 = wh.*repmat(params.nor(:,3),1,size(X2,1))./(sqrt(den+(den==0)).^3);
     NdotR = d1.*N1+d2.*N2+d3.*N3;    
     
     % diagonal part 
     A1 = (CI==CJ).*(a*(den==0)+NdotR);
     
     % w*(n'r)r_i*r_j/r^5
     A2 = 3*NdotR./(den+(den==0));
     A2(I1) = A2(I1).*d1(I1); A2(J1) = A2(J1).*d1(J1);
     A2(I2) = A2(I2).*d2(I2); A2(J2) = A2(J2).*d2(J2);
     A2(I3) = A2(I3).*d3(I3); A2(J3) = A2(J3).*d3(J3);
     
     % w*(dij*(n'r) - nirj-njri)/r^3   
     Jd = J1.*d1+J2.*d2+J3.*d3; 
     A1(I1) = A1(I1)-N1(I1).*Jd(I1); 
     A1(I2) = A1(I2)-N2(I2).*Jd(I2);
     A1(I3) = A1(I3)-N3(I3).*Jd(I3);
     
     Id = I1.*d1 + I2.*d2 + I3.*d3;
     A1(J1) = A1(J1)-N1(J1).*Id(J1);
     A1(J2) = A1(J2)-N2(J2).*Id(J2);
     A1(J3) = A1(J3)-N3(J3).*Id(J3);
     
     A = A1+A2;      
     case 'TSL_Stk_3D'    

     % create d matrices for now
     d1 = X1(:,1) - X2(:,1).';
     d2 = X1(:,2) - X2(:,2).';
     d3 = X1(:,3) - X2(:,3).';
     % Traction Stokes 
     wh = (-3/(4*pi))*wh; 
     a = params.a; 
     % Row and column tensor coords
     ci = params.ci; cj = params.cj; 
     [CJ,CI] = meshgrid(cj,ci);     
     I1 = CI==1; I2 = CI==2; I3 = CI==3;     
     J1 = CJ==1; J2 = CJ==2; J3 = CJ==3; 
     
     % N(x) dot R
     N1 = (wh./(sqrt(den+(den==0)).^3)).*params.nor(:,1);
     N2 = (wh./(sqrt(den+(den==0)).^3)).*params.nor(:,2);
     N3 = (wh./(sqrt(den+(den==0)).^3)).*params.nor(:,3);
     NdotR = (X1(:,1) - X2(:,1).').*N1+ (X1(:,2) - X2(:,2).').*N2+ (X1(:,3) - X2(:,3).').*N3;    
     % diagonal part 
     A1 = (CI==CJ).*(a*(den==0));
     
     % w*(n'r)r_i*r_j/r^5
     A2 = NdotR./(den+(den==0));
     A2(I1) = A2(I1).*d1(I1); A2(J1) = A2(J1).*d1(J1);
     A2(I2) = A2(I2).*d2(I2); A2(J2) = A2(J2).*d2(J2);
     A2(I3) = A2(I3).*d3(I3); A2(J3) = A2(J3).*d3(J3);
     
     A = A1+A2;
     case 'PSL_Stk_3D'
     % Pressure for Stokeslet (SL) 
     wh = (1/(4*pi))*wh;
     % Row and column tensor coords
     ci = params.ci; cj = params.cj; 
     [CJ,~] = meshgrid(cj,ci);       
     J1 = CJ==1; J2 = CJ==2; J3 = CJ==3; 
     
     % w*r_j/r^3
     A = wh./(sqrt(den+(den==0)).^3);
     A(J1) = A(J1).*d1(J1);
     A(J2) = A(J2).*d2(J2);
     A(J3) = A(J3).*d3(J3);
     case 'PDL_Stk_3D'
     % Pressure for Stresslet (DL) 
     wh = (1/(4*pi))*wh;
     % Row and column tensor coords
     ci = params.ci; cj = params.cj; 
     [CJ,~] = meshgrid(cj,ci);       
     J1 = CJ==1; J2 = CJ==2; J3 = CJ==3; 
     
     % N(x) dot R
     N1 = wh.*repmat(params.nor(:,1),1,size(X2,1))./(sqrt(den+(den==0)).^3);
     N2 = wh.*repmat(params.nor(:,2),1,size(X2,1))./(sqrt(den+(den==0)).^3);
     N3 = wh.*repmat(params.nor(:,3),1,size(X2,1))./(sqrt(den+(den==0)).^3);
     NdotR = d1.*N1+d2.*N2+d3.*N3; 
     
     % -w*n_j/r^3 + w*(n'r)r_i*r_j/r^5
     A = 3*NdotR./(den+(den==0));
     A(J1) = -N1(J1) + A(J1).*d1(J1);
     A(J2) = -N2(J2) + A(J2).*d2(J2);
     A(J3) = -N3(J3) + A(J3).*d3(J3); 
     case 'Rotlet_3D'
     % Rotlet Stokes 
     wh = (1/(8*pi))*wh; 
     a = params.a; 
     % Row and column tensor coords
     ci = params.ci; cj = params.cj; 
     [CJ,CI] = meshgrid(cj,ci);     
     
     % Levi-Civita index and symbol
     ijk = [0 3 2;3 0 1;2 1 0]; 
     eijk = [0 1 -1;-1 0 1;1 -1 0];
     
     K = ijk(CI + 3*(CJ-1)); 
     K = reshape(K,size(den)); 
     E = eijk(CI + 3*(CJ-1)); E = reshape(E,size(den)); 
     
     A = zeros(size(den)); 
     IK = ~(K==0);
     
     A(IK) = (E(IK).*wh(IK))./(sqrt(den(IK)+(den(IK)==0)).^3); 
     
     A(K==1) = A(K==1).*d1(K==1);
     A(K==2) = A(K==2).*d2(K==2);
     A(K==3) = A(K==3).*d3(K==3);
     
     case 'RRT_3D'
     ci = params.ci; cj = params.cj; 
     [CJ,CI] = meshgrid(cj,ci);     
     I1 = CI==1; I2 = CI==2; I3 = CI==3;     
     J1 = CJ==1; J2 = CJ==2; J3 = CJ==3; 
     A = ones(size(den)); 
     A(I1) = A(I1).*d1(I1); A(J1) = A(J1).*d1(J1);
     A(I2) = A(I2).*d2(I2); A(J2) = A(J2).*d2(J2);
     A(I3) = A(I3).*d3(I3); A(J3) = A(J3).*d3(J3);
     
     A = (den==0)+A;    
     case 'SL_K_3D'
     % Single layer elasticity 
     a = params.a; 
     ci = params.ci; cj = params.cj; 
     [CJ,CI] = meshgrid(cj,ci);     
     I1 = CI==1; I2 = CI==2; I3 = CI==3;     
     J1 = CJ==1; J2 = CJ==2; J3 = CJ==3;       
     
     % 1/r diagonal part (add dependence on elastic moduli)
     A1 = (CI==CJ).*(a*(den==0) + wh.*(1./sqrt(den + (den==0)) - (den==0)));       
     
     % r_i*r_j/r^3
     A2 = wh./(sqrt(den+(den==0)).^3);
     A2(I1) = A2(I1).*d1(I1); A2(J1) = A2(J1).*d1(J1);
     A2(I2) = A2(I2).*d2(I2); A2(J2) = A2(J2).*d2(J2);
     A2(I3) = A2(I3).*d3(I3); A2(J3) = A2(J3).*d3(J3);
     
     
     A = A1+A2; 
     %Single Layer Modified Laplace 
     case 'SL_LMOD_3D'
        % old code TESTING
        %{
        [Y_g1,  X_g1  ] = meshgrid(X2(:,1), X1(:,1));
        [Y_g2,  X_g2  ] = meshgrid(X2(:,2), X1(:,2)); 
        % Z grid
        [Y_g3,X_g3] = meshgrid(X2(:,3),X1(:,3));
        

        % den = ||X-Y||^2
        d1 = (X_g1 - Y_g1); d2 = (X_g2 - Y_g2); d3 = (X_g3 - Y_g3); 
        den_old = d1.^2 + d2.^2 + d3.^2;
        assert(isequaln(den, den_old))
        %}
        lambda=params.lambda;  
        a=params.a;
        wh = (1/4/pi)*wh;
        A = a*(den==0)+wh.*(exp(-lambda*sqrt(den))./sqrt(den + (den==0)) - (den==0));
        
     case 'DL_LMOD_3D'

        % old code TESTING
        %{
        [Y_g1,  X_g1  ] = meshgrid(X2(:,1), X1(:,1));
        [Y_g2,  X_g2  ] = meshgrid(X2(:,2), X1(:,2)); 
        % Z grid
        [Y_g3,X_g3] = meshgrid(X2(:,3),X1(:,3));
        

        % den = ||X-Y||^2
        d1 = (X_g1 - Y_g1); d2 = (X_g2 - Y_g2); d3 = (X_g3 - Y_g3); 
        den_old = d1.^2 + d2.^2 + d3.^2;
        assert(isequaln(den, den_old))
        %}
        lambda=params.lambda;
        a = params.a;
        wh = (1/(4*pi))*wh;

        % old code 
        %{
        N1 = repmat(params.nor(:,1).',size(X1,1),1);
        N2 = repmat(params.nor(:,2).',size(X1,1),1);
        N3 = repmat(params.nor(:,3).',size(X1,1),1);
        NdotR = d1.*N1+d2.*N2+d3.*N3;
        A_old = a*(den==0) + (wh.*NdotR.*exp(-lambda*sqrt(den))./sqrt(den)).*(1./den + lambda./sqrt(den)) - (den==0);
        %}
        

        % new code
        N = params.nor';
        NdotR = (X1(:,1) - X2(:,1).').*N(1,:) + (X1(:,2) - X2(:,2).').*N(2,:) + (X1(:,3) - X2(:,3).').*N(3,:);
        A = a*(den==0) + (wh.*NdotR.*exp(-lambda*sqrt(den))./sqrt(den)).*(1./den + lambda./sqrt(den)) - (den==0);


        %assert(isequaln(A_old, A))

     case 'dSL_LMOD_3D'

        % old code TESTING% old code TESTING
        %{
        [Y_g1,  X_g1  ] = meshgrid(X2(:,1), X1(:,1));
        [Y_g2,  X_g2  ] = meshgrid(X2(:,2), X1(:,2)); 
        % Z grid
        [Y_g3,X_g3] = meshgrid(X2(:,3),X1(:,3));
        

        % den = ||X-Y||^2
        d1 = (X_g1 - Y_g1); d2 = (X_g2 - Y_g2); d3 = (X_g3 - Y_g3); 
        den_old = d1.^2 + d2.^2 + d3.^2;
        assert(isequaln(den, den_old))
        %}
        lambda=params.lambda;
        a = params.a;
        wh = -(1/(4*pi))*wh;
        
        %old code
        %{
        N1 = repmat(params.nor(:,1),1,size(X2,1));
        N2 = repmat(params.nor(:,2),1,size(X2,1));
        N3 = repmat(params.nor(:,3),1,size(X2,1));
        NdotR = d1.*N1+d2.*N2+d3.*N3;    
            
        A_old = a*(den==0) + (wh.*NdotR.*exp(-lambda*sqrt(den))./sqrt(den)).*(1./den + lambda./sqrt(den)) - (den==0);
        %}

        N = params.nor;
        NdotR = (X1(:,1) - X2(:,1).').*N(:, 1) + (X1(:,2) - X2(:,2).').*N(:, 2) + (X1(:,3) - X2(:,3).').*N(:, 3);
        A = a*(den==0) + (wh.*NdotR.*exp(-lambda*sqrt(den))./sqrt(den)).*(1./den + lambda./sqrt(den)) - (den==0);

        %assert(isequaln(A_old, A))

        case 'dDL_LMOD_3D'
        lambda=params.lambda;
        a=params.a;
        wh=-(1/(4*pi))*wh;

        % old code TESTING
        %{
        [Y_g1,  X_g1  ] = meshgrid(X2(:,1), X1(:,1));
        [Y_g2,  X_g2  ] = meshgrid(X2(:,2), X1(:,2)); 
        % Z grid
        [Y_g3,X_g3] = meshgrid(X2(:,3),X1(:,3));
        

        % den = ||X-Y||^2
        d1 = (X_g1 - Y_g1); d2 = (X_g2 - Y_g2); d3 = (X_g3 - Y_g3); 
        den_old = d1.^2 + d2.^2 + d3.^2;
        assert(isequaln(den, den_old))
        
        %old code
        %computes the dot products <r,n_source>,<r,n_target>
        
        SourceN1 = repmat(params.nor(:,1).',size(X1,1),1);
        SourceN2 = repmat(params.nor(:,2).',size(X1,1),1);
        SourceN3 = repmat(params.nor(:,3).',size(X1,1),1);
        NdotRSource_old = d1.*SourceN1+d2.*SourceN2+d3.*SourceN3;
        TargN1 = repmat(params.targnor(:,1),1,size(X2,1));
        TargN2 = repmat(params.targnor(:,2),1,size(X2,1));
        TargN3 = repmat(params.targnor(:,3),1,size(X2,1));
        NdotRTarg_old = d1.*TargN1+d2.*TargN2+d3.*TargN3;  
        dotnorm_old=TargN1.*SourceN1+TargN2.*SourceN2+SourceN3.*TargN3;
        rbar_old=sqrt(den); %%%
        expy_old=exp(-lambda*rbar_old);  %%%
        A_old=(lambda^2./rbar_old + 2*lambda./rbar_old.^2 + 2./rbar_old.^3);  %%%
        B_old=(lambda./rbar_old + 1./rbar_old.^2);  %%%
        delsquared_old=-1./rbar_old.^3.*NdotRTarg_old.*NdotRSource_old + 1./rbar_old.*dotnorm_old; %%%
        A1_old =(expy_old./rbar_old.^2).*(A_old.*NdotRTarg_old.*NdotRSource_old);
        A2_old=expy_old.*B_old.*delsquared_old;
        A_old= a*(den==0) + wh.*(A1_old - A2_old) -(den==0);
        %}
        
        N = params.nor';
        T = params.targnor;
        NdotRSource = (X1(:,1) - X2(:,1).').*N(1,:)+(X1(:,2) - X2(:,2).').*N(2,:)+(X1(:,3) - X2(:,3).').*N(3,:);
        NdotRTarg = ((X1(:,1) - X2(:,1).').*T(:,1)+(X1(:,2) - X2(:,2).').*T(:,2)+(X1(:,3) - X2(:,3).').*T(:,3));
        dotnorm=T(:,1).*N(1,:)+T(:,2).*N(2,:)+T(:,3).*N(3,:);
        rbar=sqrt(den); %%%
        expy=exp(-lambda*rbar);  %%%
        A=(lambda^2./rbar + 2*lambda./rbar.^2 + 2./rbar.^3);  %%%
        B=(lambda./rbar + 1./rbar.^2);  %%%
        delsquared=-1./rbar.^3.*NdotRTarg.*NdotRSource + 1./rbar.*dotnorm; %%%
        A1 =(expy./rbar.^2).*(A.*NdotRTarg.*NdotRSource);
        A2=expy.*B.*delsquared;
        A= a*(den==0) + wh.*(A1 - A2) -(den==0);
        
        %{
        assert(isequaln(rbar, rbar_old));
        assert(isequaln(expy, expy_old));
        assert(isequaln(delsquared, delsquared_old));
        assert(isequaln(B, B_old));
        assert(isequaln(A1, A1_old));
        assert(isequaln(A2, A2_old));
        assert(isequaln(A, A_old))
        %}

     case 'Distance'
        A = sqrt(den); 
     case 'zeros'
         A = zeros(size(den)); 
     case 'Fun_3D'
        fun =params.fun; 
        A  =wh.*fun(den);  
    case 'Kernel_Name'
        % User-defined Kernel template K(x,y) = I + b(X)K(|X-Y|)c(Y) 
        A = (den==0) + 0.5*wh.*b(X_g1,X_g2).*K(den).*c(Y_g1,Y_g2);    
     end
 end
 
end
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 % b function for NTI Laplace example
function b = b_func(xx,x0,alpha)

if alpha==0
    b = ones(size(xx{1})); 
else
    dd2 = ones(size(xx{1})); 
    for i=1:length(xx)
        dd2 = dd2 + (xx{i}-x0(i)).*(xx{i}-x0(i)); 
    end
    
    b   = 1 + alpha*exp( - dd2);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [b,c] = LOCAL_get_bc(Xg,Yg,params,opts)

if params.sym == 0 && params.kh>0
    % Option for non-symmetric Lippmann-Schwinger 
    if params.proxy <= 0
    b = Bump_fun(Xg,opts); 
    else
    b = ones(size(Xg{1})); 
    end
    
    if params.proxy >= 0 
        c = ones(size(Yg{1})); 
        %c = Bump_function(Y_g1,Y_g2,opts);
    else
        c = ones(size(Yg{1})); 
    end
elseif params.transinv == 0 
    % Option for symmetric Laplace / Lippmann-Schwinger 
    if params.proxy <= 0
        if params.kh>0
            b = sqrt(Bump_fun(Xg,opts)); 
        else
            b = b_func(Xg,[0.3 0.6 -0.4],0.25); 
        end
    else
        b = ones(size(Xg{1})); 
    end
    
    if params.proxy >=0
        if params.kh>0
            c = sqrt(Bump_fun(Yg,opts));
        else
            c = b_func(Yg,[0.3 0.6 -0.4],0.25);      
        end
    else
        c = ones(size(Yg{1})); 
    end  
else
    % Else, b and c are ones
    b = ones(size(Xg{1})); 
    c = ones(size(Yg{1})); 
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [params,X1,X2] = set_params(params,X1,X2)
if strcmp(params.keval,'ind')
   II = X1; JJ = X2; 
   X1 = params.X(II,:); X2 = params.X(JJ,:); 
   if isfield(params,'nor')
       if strcmp(params.flag_pot,'dSL_Stk_3D') || strcmp(params.flag_pot,'dSL_L_3D') || strcmp(params.flag_pot,'TSL_Stk_3D')
           params.nor = params.nor(II,:); 
       else
           params.nor = params.nor(JJ,:); 
       end
   end
   
   if isfield(params,'W2')
      params.W2 = params.W2(JJ);  
   end
   
   if isfield(params,'ci')
       params.ci = params.ci(II); 
       params.cj = params.cj(JJ);    
   end
end
end
