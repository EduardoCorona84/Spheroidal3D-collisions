function [CM, CM_trans, CM_rot] = RBM_completion(pars, Y, ns)
    %{
        Calculates the completion term that corresponds to rigid body
        motion.

        CM*[sigma_x ; sigma_y; sigma_z].
    %}
    p = pars.p; np = 2*p*(p+1);

    % Unit quadrature weights
    [~, gwt]=g_grid(p+1);
    wt = pi/p*repmat(gwt', 2*p, 1)./sin(gl_grid(p));
    wt = wt(:);

    CM_cells = cell(1, ns);
    for i=1:ns
        if ~pars.oblate(i)
            Yi = prolate_spheroid_shape(p,pars.u0(i),pars.a(i),'cart');
        else
            Yi = oblate_spheroid_shape(p,pars.u0(i),pars.a(i),'cart');
        end
        ci = pars.centers(i,:);

        % Surface area weights
        Sns = SurfaceSph(Yi(:));
        W = Sns.geoProp.W;

        % Scaled quadrature weights
        Wns = W.*wt;

        % Calculate moment of inertia tensor \tau
        R = Yi - ci;
        M = R' * (Wns .* R);
        tau = trace(M) * eye(3) - M;
        tauTest = [
            Wns.'*(R(:,2).^2 + R(:,3).^2) -1*Wns.'*(R(:,1).*R(:,2)) -1*Wns.'*(R(:,1).*R(:,3));
            -1*Wns.'*(R(:,1).*R(:,2)) Wns.'*(R(:,1).^2 + R(:,3).^2) -1*Wns.'*(R(:,2).*R(:,3));
            -1*Wns.'*(R(:,1).*R(:,3)) -1*Wns.'*(R(:,2).*R(:,3)) Wns.'*(R(:,1).^2 + R(:,2).^2);
        ];

        % Rotational portion
        rotational_integral_term = [
            zeros(1, np),       -1*(Wns.*R(:,3)).', (Wns.*R(:,2)).';
            (Wns .* R(:,3)).',    zeros(1, np),     -1*(Wns.*R(:,1)).';
            -1*(Wns.*R(:,2)).',   (Wns.*R(:,1)).',    zeros(1, np)
        ]; % 3 x 3np (input is density)

        % This represents the cross-product with the integral (accounting 
        % for the anti-commutativity),
        cross_prod_term = -1*[
            zeros(np, 1), -1*R(:,3),    R(:,2);
            R(:,3),       zeros(np, 1), -1*R(:,1);
            -1*R(:,2),    R(:,1),       zeros(np, 1)
        ]; % 3np x 3

        CM_rot = cross_prod_term * (tau \ rotational_integral_term); % 3np x 3np

        % Calculate translational portion
        surf_area = sum(Wns);
        T = (1/surf_area) * ones(np, 1) * Wns';
        CM_trans = kron(eye(3), T); % 3np x 3np

        CM_cells{i} = CM_trans + CM_rot;
    end
    CM = blkdiag(CM_cells{:});
end

  