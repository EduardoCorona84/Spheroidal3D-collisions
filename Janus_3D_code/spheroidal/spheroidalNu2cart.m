function [Nu,recons] = spheroidalNu2cart(nu_sph,strg,a,oblate)
    % nu_sph is the nu in spheroidal coordinates
    % strg are the spheroidal coordinates of the target in the frame of
    % source.
    % a is the eccentricity of the source sphere.
    if nargin<1
        test2();
        return;
    end
    if nargin==3
        oblate=0;
    end
    if size(nu_sph,1)==3
        nu = nu_sph';
    elseif size(nu_sph,2)==3
        nu=nu_sph;
    else
        fprintf("\n inputting Nu should have dimension 3.")
    end
    un=strg(:,1); vn=strg(:,2); phin=strg(:,3);
    if ~oblate
        h_u = a.*sqrt((un.^2-vn.^2)./(un.^2-1));
        h_v = a.*sqrt((un.^2-vn.^2)./(1-vn.^2));
        h_phi = a.*sqrt((un.^2-1).*(1-vn.^2));
        eu=[a.*un.*sqrt((1-vn.^2)./(un.^2-1)).*cos(phin),a.*un.*sqrt((1-vn.^2)./(un.^2-1)).*sin(phin),a.*vn];
        ev=[-a.*vn.*sqrt((un.^2-1)./(1-vn.^2)).*cos(phin),-a.*vn.*sqrt((un.^2-1)./(1-vn.^2)).*sin(phin),a.*un];
        ephi=[-a.*sqrt((un.^2-1).*(1-vn.^2)).*sin(phin),a.*sqrt((un.^2-1).*(1-vn.^2)).*cos(phin),zeros(size(un))];
    else
        h_u = a.*sqrt((un.^2+vn.^2)./(un.^2+1));
        h_v = a.*sqrt((un.^2+vn.^2)./(1-vn.^2));
        h_phi = a.*sqrt((un.^2+1).*(1-vn.^2));
        eu=[a.*un.*sqrt((1-vn.^2)./(un.^2+1)).*cos(phin),a.*un.*sqrt((1-vn.^2)./(un.^2+1)).*sin(phin),a.*vn];
        ev=[-a.*vn.*sqrt((un.^2+1)./(1-vn.^2)).*cos(phin),-a.*vn.*sqrt((un.^2+1)./(1-vn.^2)).*sin(phin),a.*un];
        ephi=[-a.*sqrt((un.^2+1).*(1-vn.^2)).*sin(phin),a.*sqrt((un.^2+1).*(1-vn.^2)).*cos(phin),zeros(size(un))];
    end
    eu = eu./h_u;
    ev = ev./h_v;
    ephi = ephi./h_phi;
    if size(nu_sph)~=size(strg)
        fprintf("\n Size of nu not equal to size of Strg.\n");
    end
    Nu = nu(:,1).*eu+nu(:,2).*ev+nu(:,3).*ephi;
    recons = [diag(Nu*eu'), diag(Nu*ev'), diag(Nu*ephi')];

    function test()
        nu_sph_test=[1,0,0;0,1,0;0,0,1];
        u0_test=2/sqrt(3); a_test=1/u0_test;
        u_r = [0.5;0.6;0.7]+u0_test; v_r = repmat(0.58,3,1); phi_r = repmat(pi/7,3,1);
        Xtrg_test = [a_test.*sqrt(u_r.^2-1).*sqrt(1-v_r.^2).*cos(phi_r),a_test.*sqrt(u_r.^2-1).*sqrt(1-v_r.^2).*sin(phi_r),a_test.*u_r.*v_r]; 
        Xpos = [a_test.*sqrt(u0_test.^2-1).*sqrt(1-v_r.^2).*cos(phi_r),a_test.*sqrt(u0_test.^2-1).*sqrt(1-v_r.^2).*sin(phi_r),a_test.*u0_test.*v_r];
        display(Xpos)
        Spt_test = cart2spheroidal(Xtrg_test,a_test);
        [Nu_test,recons_test] = spheroidalNu2cart(nu_sph_test,Spt_test,a_test);
        display(Nu_test);
        display(recons_test);
        
        Ssp = prolate_spheroid_shape(8,u0_test,a_test);
        plot3(Ssp(:,1), Ssp(:,2),Ssp(:,3)); hold on;
        quiver3(Xpos(1,1),Xpos(1,2),Xpos(1,3), Nu_test(1,1),Nu_test(1,2),Nu_test(1,3),'r'); 
        quiver3(Xpos(2,1),Xpos(2,2),Xpos(2,3), Nu_test(2,1),Nu_test(2,2),Nu_test(2,3),'r'); 
        quiver3(Xpos(3,1),Xpos(3,2),Xpos(3,3), Nu_test(3,1),Nu_test(3,2),Nu_test(3,3),'r'); hold off; axis equal
    end

    function test2()
        u0_test=2/sqrt(3); a_test=1/u0_test;
        Xtrg=prolate_spheroid_shape(8,u0_test,a_test);
        Strg=cart2spheroidal(Xtrg,a_test);
        % Xtrg_dist=Xtrg+2.*ones(size(Xtrg));
        % Strg_dist=cart2spheroidal(Xtrg_dist,a_test);
        Shape_Strg=SurfaceSph(Xtrg);
        nu_cart = reshape(Shape_Strg.geoProp.nor.to_array,[],3);
        quiver3(Xtrg(:,1),Xtrg(:,2),Xtrg(:,3),nu_cart(:,1),nu_cart(:,2),nu_cart(:,3),'r');

        nu_sph_test=repmat([1,0,0],size(Xtrg,1),1);
        [Nu_test,recons_test]=spheroidalNu2cart(nu_sph_test,Strg,a_test);
        figure;
        quiver3(Xtrg(:,1),Xtrg(:,2),Xtrg(:,3),Nu_test(:,1),Nu_test(:,2),Nu_test(:,3),'r');
        display(norm(abs(nu_cart-Nu_test)));
    end
end