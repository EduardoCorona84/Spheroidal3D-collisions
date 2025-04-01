%{
    Test code for Laplace to Stokes single layer with two prolate
    spheroids.

    The density, sigma, is randomized.
%}

clear;

DESIRED_TOL = 1e-6;
CHECK_LAPLACE_FLAG = false;
CHECK_STOKES_FLAG = true;
CHECK_RUN_OLD_CODE = false;

params=SpheroidalParameters; 
params.matvec_eta=10; 
params.u0=1.1;
params.a=1/1.1;
params.oblate=0; 
params.centers=[0 0 0];
p=8; np=2*p*(p+1);

sigma_x=rand(np,1)+0.5;
sigma_y=rand(np,1)-0.5;
sigma_z=rand(np,1);

% Get Cartesian coordinates from the two spheroids
Xtrg=prolate_spheroid_shape(p,1.2,1/1.2)+[3,3,1];
Xself=prolate_spheroid_shape(p, params.u0, params.a);

Sns = SurfaceSph(Xself);
integrate = @(f) integrateOverS(Sns, f);

Ntrg = size(Xtrg, 1);
Nself = size(Xself, 1);

% Pre-calculate Rvec and r for all target-source pairs
Rvec_x = Xtrg(:, 1)' - Xself(:, 1);
Rvec_y = Xtrg(:, 2)' - Xself(:, 2);
Rvec_z = Xtrg(:, 3)' - Xself(:, 3);
r = sqrt(Rvec_x.^2 + Rvec_y.^2 + Rvec_z.^2);

if CHECK_LAPLACE_FLAG
    sigma_list = [sigma_x, sigma_y, sigma_z]; % Size Nself x 3
    Imat = eye(3);
    
    fprintf('\n--- Checking Laplace potentials... ---\n');
    for sigma_ind = 1:3
        sig = sigma_list(:, sigma_ind);
        params.sigma = sig;
        params.get_shc;
    
        % Create kernel of Laplace SLP
        Sl_sig_vec(1, :, sigma_ind) = (1/4/pi) * integrate(sig ./ r);
    
        % Calculate MatVec for all target points
        SL_matvec_vec(:, 1, sigma_ind) = spheroidalMatVec(params, 'SL', Xtrg);
    
        fprintf("SL sigma_%d vs matvec max abs error: %e\n", sigma_ind, ...
                max(norm( Sl_sig_vec(1,:,sigma_ind)' - SL_matvec_vec(:, 1, sigma_ind) )) );
    
        for j = 1:3 % Derivative direction index (dx, dy, dz)
            % Select correct Rvec component based on j
            if j == 1; Rvec_comp = Rvec_x;
            elseif j == 2; Rvec_comp = Rvec_y;
            else; Rvec_comp = Rvec_z;
            end
    
            dSl_sig_vec(1, :, sigma_ind, j) = -1/(4*pi) * integrate(sig .* Rvec_comp ./ (r.^3));
    
            nu_stacked_vec = repmat(Imat(j,:), size(Xtrg,1), 1); % Create a normal vector for every target point
            dSL_matvec_vec(:, 1, sigma_ind, j) = spheroidalMatVec(params, 'SP', Xtrg, nu_stacked_vec);
    
            fprintf("SP sigma_%d dx%d vs matvec max abs error: %e\n", sigma_ind, j, ...
                    max(norm( dSl_sig_vec(1, :, sigma_ind, j)' - dSL_matvec_vec(:, 1, sigma_ind, j) )) );
        end
    end
    fprintf('------ Done... ------\n');
end

if CHECK_STOKES_FLAG
    % Precomputations
    r_inv = 1 ./ r;
    r_inv_3 = r_inv .^ 3;

    % Integrand computations
    Fx      = sigma_x .* r_inv;
    F11x    = sigma_x .* Rvec_x .* Rvec_x .* r_inv_3;
    F12y    = sigma_y .* Rvec_x .* Rvec_y .* r_inv_3;
    F13z    = sigma_z .* Rvec_x .* Rvec_z .* r_inv_3;
    integrand_SLx = Fx + F11x + F12y + F13z;

    Fy      = sigma_y .* r_inv;
    F21x    = sigma_x .* Rvec_y .* Rvec_x .* r_inv_3;
    F22y    = sigma_y .* Rvec_y .* Rvec_y .* r_inv_3;
    F23z    = sigma_z .* Rvec_y .* Rvec_z .* r_inv_3;
    integrand_SLy = Fy + F21x + F22y + F23z;

    Fz      = sigma_z .* r_inv;
    F31x    = sigma_x .* Rvec_z .* Rvec_x .* r_inv_3;
    F32y    = sigma_y .* Rvec_z .* Rvec_y .* r_inv_3;
    F33z    = sigma_z .* Rvec_z .* Rvec_z .* r_inv_3;
    integrand_SLz = Fz + F31x + F32y + F33z;

    % Integrate
    SLx = (1/(8*pi) * integrate(integrand_SLx))';
    SLy = (1/(8*pi) * integrate(integrand_SLy))';
    SLz = (1/(8*pi) * integrate(integrand_SLz))';

    % Now, actually calculate result from L2Stk and compare.
    target_pts = cell(1, 1);
    target_pts{1} = Xtrg;
    [L2Stkx, L2Stky, L2Stkz] = L2Stk(target_pts, params, sigma_x, sigma_y, sigma_z, 1);
    fprintf("\n inf error Stokes in x: %e, in y: %e, in z: %e\n", max(norm(L2Stkx{1}-SLx)), max(norm(L2Stky{1}-SLy)),max(norm(L2Stkz{1}-SLz)));
end


%%%%%%%%%%%%
% OLD CODE %
%%%%%%%%%%%%
if CHECK_RUN_OLD_CODE
    for trgind=1:size(Xtrg,1)
        % was having trouble with Rvec for multiple target points. So for
        % testing, loop over different test points. However, no restrictions on
        % test points for actual L2Stk function, so can vectorize.
        
        Xeval=Xtrg(trgind,:);
        Rvec=[Xeval(1)-Xself(:,1) Xeval(2)-Xself(:,2) Xeval(3)-Xself(:,3)];
        r=(vecnorm(Rvec'))';
    
        %% first check Laplace
        %{
        sigma_list=[sigma_x sigma_y sigma_z];
        Imat=eye(3);
        for sigma_ind=1:3
            sig=sigma_list(:,sigma_ind);
            params.sigma=sig;
            params.get_shc;
        
            SL_sig=1/4/pi.*integrateOverS(Sns,sig./r);
            SL_matvec=spheroidalMatVec(params,'SL',Xeval);
            fprintf("\n SL sigma_i vs matvec inf error %e", max(abs(SL_sig-SL_matvec)));
        
            for j=1:3
                dSl_sig=-1/4/pi.*integrateOverS(Sns,sig.*Rvec(:,j)./(r.^3));
                dSL_matvec=spheroidalMatVec(params,'SP',Xeval,Imat(j,:));
                fprintf("\n SP sigma_i dxj inf error %e",max(abs(dSl_sig-dSL_matvec)));
            end
        end
        %}
    
        %% Next check Stokes
        % %{
        % i=1
        F1=sigma_x./r;
        
        F31x=sigma_x.*Rvec(:,1).*Rvec(:,1)./(r.^3);
        F32y=sigma_y.*Rvec(:,1).*Rvec(:,2)./(r.^3);
        F13z=sigma_z.*Rvec(:,1).*Rvec(:,3)./(r.^3);
        
        SLx=1/8/pi.*integrateOverS(Sns,F1+F31x+F32y+F13z);
        
        % i=2
        F1=sigma_y./r;
        
        F12x=sigma_x.*Rvec(:,2).*Rvec(:,1)./(r.^3);
        F32y=sigma_y.*Rvec(:,2).*Rvec(:,2)./(r.^3);
        F33z=sigma_z.*Rvec(:,2).*Rvec(:,3)./(r.^3);
        
        Sns=SurfaceSph(Xself);
        SLy=1/8/pi.*integrateOverS(Sns,F1+F12x+F32y+F33z);
        
        % i=3
        F1=sigma_z./r;
        
        F13x=sigma_x.*Rvec(:,3).*Rvec(:,1)./(r.^3);
        F23y=sigma_y.*Rvec(:,3).*Rvec(:,2)./(r.^3);
        F33z=sigma_z.*Rvec(:,3).*Rvec(:,3)./(r.^3);
        
        Sns=SurfaceSph(Xself);
        SLz=1/8/pi.*integrateOverS(Sns,F1+F13x+F23y+F33z);
        
        % calculated from laplace operators.
        [theta,phi]=gl_grid(p);
        % v=cos(theta);
        target_pts = cell(1, 1);
        target_pts{1} = Xeval;
        [L2Stkx,L2Stky,L2Stkz]=L2Stk(target_pts,params,sigma_x,sigma_y,sigma_z,1);
        fprintf("\n inf error Stokes in x: %e, in y: %e, in z: %e\n", max(abs(L2Stkx{1}-SLx)),max(abs(L2Stky{1}-SLy)),max(abs(L2Stkz{1}-SLz)));
        % %}
    end
end