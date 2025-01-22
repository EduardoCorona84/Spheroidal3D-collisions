function A = Kernel_Eval_vec(X1,X2,params)
%
% Modified file from HSSDirectSolver. 
% Copyright (C) 2013 Eduardo Corona, Denis Zorin, Per Gunnar Martinsson
% See <COPYRIGHT_NOTICE.txt> for more details.
%
%
%     FUNCTION CALL:
%         A = Kernel_Eval(X1,X2,params)
%
%     top - for nonsym multiplication by b and c
%     DESCRIPTION:
%         This function takes 2D point arrays X1 and X2, and returns the 
%         vector with entries A(i) = K[X1(i,:),X2(i,:)].  
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
    wh = params.W2.';        
end

% Determine functions b(x) and c(y)
%[b,c] = LOCAL_get_bc(X1,X2,params,4);

% 2D (plane or curves on plane)
 if params.dim == 2
     % den = ||X-Y||^2
     den = (X1(:,1) - X2(:,1)).^2 + (X1(:,2) - X2(:,2)).^2;   
    
    switch params.flag_pot
     % 2D Single Layer Helmholtz
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
    % 3D Single Layer Helmholtz    
    case 'SL_H_3D'
        kh = params.kh; 
        C = exp(1i*kh);
        A = (den==0) + wh*kh*(exp(1i*kh*sqrt(den + (den==0)))./sqrt(den + (den==0)) - C*(den==0)); 
    % 2D Single Layer Laplace    
    case 'SL_L_2D'
        wh = (1/2/pi)*wh; 
        A = (den==0) + 0.5*wh.*b.*log((den + (den==0))).*c;
    % 3D Single Layer Laplace    
    case 'SL_L_3D'
        wh = (1/4/pi)*wh; 
        A = (den==0) + wh.*(1./sqrt(den + (den==0)) - (den==0)); 
    case 'DO_L_2D_5' 
        den = (1/h)*sqrt(den);     
        A = -4*(den==0)+1*(den==1); 
    case 'DO_L_2D_9'
        den = (1/h)*sqrt(den);      
        A = -(20/6)*(den==0)+(4/6)*(den==1)+(1/6)*(den==sqrt(2));       
    % 2D Single Layer Yukawa    
    case 'SL_Y_2D'
        C = besselk(0,params.kh);  
        A = (den==0) + wh*(0.5/pi)*b.*(besselk(0,params.kh*(den + (den==0))) - C*(den==0)).*c; 
    %3D Single Layer Yukawa    
    case 'SL_Y_3D'
        A = (den==0) + wh*b.*(exp(-params.kh*den)./sqrt(den + (den==0)) - (den==0)).*c;  
    case 'Fun_2D'
        fun =params.fun; 
        A  =wh.*fun(den);  
    case 'Kernel_Name'
        % User-defined Kernel template K(x,y) = I + b(X)K(|X-Y|)c(Y) 
        A = (den==0) + 0.5*wh.*b(X_g1,X_g2).*K(den).*c(Y_g1,Y_g2);
    end
 % 3D (surfaces and volume)   
 else
     
     % den = ||X-Y||^2
     den = (X1(:,1) - X2(:,1)).^2 + (X1(:,2) - X2(:,2)).^2 + (X1(:,3) - X2(:,3)).^2;   
    
    switch params.flag_pot 
    case 'SL_L_3D'
        wh = (1/(4*pi))*wh;
        A = (den==0) + wh.*b.*(1./sqrt(den + (den==0)) - (den==0)).*c;
    case 'DL_L_3D'   
        % Double layer Laplace
        a = params.a; 
        wh = (1/(4*pi))*wh;
        NdotR = (X1(:,1) - X2(:,1)).*params.nor(:,1)+(X1(:,2) - ...
            X2(:,2)).*params.nor(:,2)+(X1(:,3) - X2(:,3)).*params.nor(:,3);           
        A = a*(den==0) + wh.*NdotR.*(1./(den+(den==0)).^(3/2) - (den==0)); 
    case 'dSL_L_3D'   
        % Double layer Laplace
        a = params.a;
        wh = (-1/(4*pi))*wh;
        NdotR = (X1(:,1) - X2(:,1)).*params.nor(:,1)+(X1(:,2) - ...
            X2(:,2)).*params.nor(:,2)+(X1(:,3) - X2(:,3)).*params.nor(:,3);           
        A = a*(den==0) + wh.*NdotR.*(1./(den+(den==0)).^(3/2) - (den==0));     
     % 3D Helmholtz    
    case 'SL_H_3D'
        kh = params.kh; 
        C = exp(1i*kh);
        A = (den==0) + kh*wh.*b.*(exp(1i*kh*sqrt(den + (den==0)))./sqrt(den + (den==0)) - C*(den==0)).*c; 
    % Differential Operator Laplace (7pt stencil)    
    case 'DO_L_3D_7'   
        den = (1/h)*sqrt(den);     
        A = -6*(den==0)+1*(den==1);    
     case 'SL_Stk_3D' 
     % Single layer Stokes 
        wh = (1/(8*pi*params.mu))*wh; 
        a = params.a; 
        % Row and column tensor coords
        ci = params.ci; cj = params.cj; 
        
        % Diagonal part (~1/r)
        A = (ci==cj).*(a*(den==0) + wh.*(1./sqrt(den + (den==0)) - (den==0))); 
        
        m = size(X1,1); 
        indi = m*(ci-1)+(1:m).'; indj = m*(cj-1)+(1:m).'; 
     
        A = A + wh.*((X1(indi)-X2(indi)).*(X1(indj)-X2(indj)))./(sqrt(den+(den==0)).^3); 
     case 'RPY_Stk_3D'
     % Rotne-Prager-Yamakawa tensor (assumes a=1 for now)
     
     wh = (1/(pi*params.mu))*wh; 
     a = params.a; 
     % Row and column tensor coords
     ci = params.ci; cj = params.cj;
     iidx = ci==cj; 
     m = size(X1,1); 
     indi = m*(ci-1)+(1:m).'; indj = m*(cj-1)+(1:m).';
     rxr = ((X1(indi)-X2(indi)).*(X1(indj)-X2(indj)));
     % Create index to separate den<2 from den>=2
     den = sqrt(den);
     nr = den<2; fr=~nr;
     
     A = zeros(size(den));  
     
     % Anear = (k*T/6*pi*mu*a)[(1-(9/32)*(r/a)) I + (3/32a)(rxr/r)]
     dnr = den(nr); 
     A(nr) = iidx(nr).*(a*(dnr==0) + (1/6)*wh(nr).*(1-(9/32)*dnr)) ...
         + (1/64)*wh(nr).*rxr(nr)./((dnr+(dnr==0)));
     
     % Afar = (k*T/8*pi*mu)[(1/r + 2a^2/3r^3) I + (rxr/r.^3 - 2a^2 rxr/r^5)]
     dfr = den(fr); 
     A(fr) = (1/8)*wh(fr).*(iidx(fr).*(1./dfr + (2/3)./(dfr.^3)) ...
         + rxr(fr).*(1./dfr.^3 - 2./dfr.^5));
     
     case 'DL_Stk_3D'  
        % Double layer Stokes 
        wh = (3/(4*pi))*wh; 
        a = params.a; 
        % Row and column tensor coords
        ci = params.ci; cj = params.cj; 
     
        % N(y) dot R
        NdotR = (X1(:,1) - X2(:,1)).*params.nor(:,1)+(X1(:,2) - ...
            X2(:,2)).*params.nor(:,2)+(X1(:,3) - X2(:,3)).*params.nor(:,3); 
     
        % diagonal part 
        A = (ci==cj).*(a*(den==0));
     
        % r_i*r_j/r^5
        m = size(X1,1); 
        indi = m*(ci-1)+(1:m).'; indj = m*(cj-1)+(1:m).'; 
     
        A = A + NdotR.*wh.*((X1(indi)-X2(indi)).*(X1(indj)-X2(indj)))./(sqrt(den+(den==0)).^5);    
     case 'dSL_Stk_3D'    
     % dS/dn_x Stokes 
        wh = -(1/(8*pi))*wh;        
        a = params.a; 
        % Row and column tensor coords
        ci = params.ci; cj = params.cj; 
     
        % N(x) dot R
        NdotR = (X1(:,1) - X2(:,1)).*params.nor(:,1)+(X1(:,2) - ...
            X2(:,2)).*params.nor(:,2)+(X1(:,3) - X2(:,3)).*params.nor(:,3); 
     
        m = size(X1,1); 
        indi = m*(ci-1)+(1:m).'; indj = m*(cj-1)+(1:m).'; 
        Zi = X1(indi)-X2(indi); Zj = X1(indj)-X2(indj); 
        rho = sqrt(den+(den==0)); 
        
        %dij(1+w(n'r)/r^3) - w*(nirj + njri)/r^3) + 3*(n'r)*r_i*r_j/r^5
        A = (ci==cj).*(a*(den==0))+ wh.*(((ci==cj).*NdotR ...
            - params.nor(indi).*Zj - params.nor(indj).*Zi)./(rho.^3)+...
            3*NdotR.*Zi.*Zj./(rho.^5));
    case 'TSL_Stk_3D'    
     % Traction Stokes 
        wh = (3/(4*pi))*wh;     
        a = params.a; 
        % Row and column tensor coords
        ci = params.ci; cj = params.cj; 
     
        % N(x) dot R
        NdotR = (X1(:,1) - X2(:,1)).*params.nor(:,1)+(X1(:,2) - ...
            X2(:,2)).*params.nor(:,2)+(X1(:,3) - X2(:,3)).*params.nor(:,3); 
     
        m = size(X1,1); 
        indi = m*(ci-1)+(1:m).'; indj = m*(cj-1)+(1:m).'; 
        Zi = X1(indi)-X2(indi); Zj = X1(indj)-X2(indj); 
        rho = sqrt(den+(den==0)); 
        
        %w*(n'r)*r_i*r_j/r^5
        A = (ci==cj).*(a*(den==0))+ wh.*NdotR.*Zi.*Zj./(rho.^5); 
    case 'PSL_Stk_3D'
     % Pressure for Stokeslet (SL) 
     wh = (1/(4*pi))*wh;
     % Row and column tensor coords
     cj = params.cj; 
     m = size(X1,1); 
     indj = m*(cj-1)+(1:m).'; 
     Zj = X1(indj)-X2(indj); 
     
     % w*r_j/r^3
     A = wh.*Zj./(sqrt(den+(den==0)).^3);
     case 'PDL_Stk_3D'
     % Pressure for Stresslet (DL) 
     wh = (1/(4*pi))*wh;
     % Column tensor coords
     cj = params.cj; 
     
     % N(x) dot R
     NdotR = (X1(:,1) - X2(:,1)).*params.nor(:,1)+(X1(:,2) - ...
     X2(:,2)).*params.nor(:,2)+(X1(:,3) - X2(:,3)).*params.nor(:,3); 
     
     m = size(X1,1); 
     indj = m*(cj-1)+(1:m).'; 
     Zj = X1(indj)-X2(indj); 
     rho = sqrt(den+(den==0)); 
     
     % -w*n_j/r^3 + w*(n'r)r_i*r_j/r^5
     A = -wh.*params.nor(indj)./rho.^3 + 3*wh.*Zj.*NdotR./rho.^5;
     A(J1) = -N1(J1) + A(J1).*d1(J1);
     A(J2) = -N2(J2) + A(J2).*d2(J2);
     A(J3) = -N3(J3) + A(J3).*d3(J3);     
    % Elasticity (Kelvin) 3D Single Layer     
    case 'SL_K_3D'
        ci = params.ci; cj = params.cj; 
        
        % Diagonal part (~1/r)
        A = (ci==cj).*((den==0) + wh.*(1./sqrt(den + (den==0)) - (den==0))); 
        
        m = size(X1,1); 
        indi = m*(ci-1)+(1:m).'; indj = m*(cj-1)+(1:m).'; 
     
        A = A + wh.*((X1(indi)-X2(indi)).*(X1(indj)-X2(indj)))./(sqrt(den+(den==0)).^3);
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
function [b,c] = LOCAL_get_bc(X1,X2,params,opts)
    
h = params.h; 
dim = params.dim; 
Xg = cell(dim,1); Yg=Xg; 
for i=1:dim
    Xg{i} = X1(:,i); 
    Yg{i} = X2(:,i);   
end

if params.transinv == 1   
    %b and c are ones
    b = ones(size(Xg{1})); 
    c = ones(size(Yg{1})); 
elseif params.sym == 0 && params.kh>0
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
else
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
