function [S,recons] = cartNu2spheroidal(nu_cart,strg,a,oblate)
    if nargin<1
        test(1);
        return;
    end
    if nargin<4
        oblate=0;
    end
    if size(nu_cart,2)==3
        nu=nu_cart;
    elseif size(nu_cart,1)==3
        nu=nu_cart';
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
    if size(nu)~=size(strg)
        fprintf("\n Size of nu not equal to size of Strg.\n");
    end
    S = [diag(nu*eu'), diag(nu*ev'), diag(nu*ephi')];
    recons = S(:,1).*eu+S(:,2).*ev+S(:,3).*ephi;

    function test(obl)
        u0_test=2/sqrt(3);
        if ~obl
            a_test=1/u0_test;
            Xtrg=prolate_spheroid_shape(8,u0_test,a_test);
        else
            a_test=1/sqrt(1+u0_test^2);
            Xtrg=oblate_spheroid_shape(8,u0_test,a_test);
        end
        Strg=cart2spheroidal(Xtrg,a_test,obl);
        Shape_Strg=SurfaceSph(Xtrg);
        nu_cart_test = reshape(Shape_Strg.geoProp.nor.to_array,[],3);
        % quiver3(Xtrg(:,1),Xtrg(:,2),Xtrg(:,3),nu_cart_test(:,1),nu_cart_test(:,2),nu_cart_test(:,3),'r');
        
        nu_sph_test=repmat([1,0,0],size(Xtrg,1),1);
        [Nu_test,recons_test] = cartNu2spheroidal(nu_cart_test,Strg,a_test,obl);
        display(norm(abs(nu_sph_test-Nu_test)));
        display(norm(abs(nu_cart_test-recons_test)));
    end

    function test2()
        u0_test=2/sqrt(3); a_test=1/u0_test;
        Xtrg=prolate_spheroid_shape(8,u0_test,a_test);
        Strg=cart2spheroidal(Xtrg,a_test);
        % Xtrg_dist=Xtrg+2.*ones(size(Xtrg));
        % Strg_dist=cart2spheroidal(Xtrg_dist,a_test);
        Shape_Strg=SurfaceSph(Xtrg);
        nu_cart_test = reshape(Shape_Strg.geoProp.nor.to_array,[],3);
        quiver3(Xtrg(:,1),Xtrg(:,2),Xtrg(:,3),nu_cart_test(:,1),nu_cart_test(:,2),nu_cart_test(:,3),'r');

        nu_sph_test=repmat([1,0,0],size(Xtrg,1),1);
        [Nu_test,recons_test]=cartNu2spheroidal(nu_cart_test,Strg,a_test);
        display(norm(abs(nu_sph_test-Nu_test)));
    end
end
