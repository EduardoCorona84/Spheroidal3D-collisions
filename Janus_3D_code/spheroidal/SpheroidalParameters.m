classdef SpheroidalParameters < matlab.mixin.Copyable
    properties 
        sigma
        sigma_coefficients
        u0
        a
        oblate = false % boolean for oblate spheroid coordinates
        
        % [x,y,z] coordinates of spheroid centers
        centers=[0 0 0]

        % Convention is that, to rotate spheroid, we rotate by theta 
        % (around the y-axis) and then phi (around the z-axis).
        thetas=[0]
        phis=[0]
        
        % alternatively, store rotation matrix for each spheroid in system
        % rather than axis-wise rotations.
        Rmat=[];

        isReal = false
        matvec_eta %cutoff distance parameter for matvec

        etaS
        etaC
    end
    properties
        p = 0
    end
 
    methods
        % Automatically calculate 'p' once sigma or sigma_coefficients is
        % provided
        function p = get.p(obj)
            if isempty(obj.sigma) && isempty(obj.sigma_coefficients)
                p=0;
            elseif ~isempty(obj.sigma)
                np = size(obj.sigma,1);
                p=sqrt(np/2+1/4)-1/2;
                if p~=fix(p)
                    error("size of sigma must match discretization on a Gauss-Legendre grid.")
                end
            elseif isempty(obj.sigma) && ~isempty(obj.sigma_coefficients)
                sp=size(obj.sigma_coefficients,1);
                p=sqrt(sp)-1;
                if p~=fix(p)
                    error("size of sigma must match discretization on a Gauss-Legendre grid.")
                end
            end
            
        end

        %Calculate spherical/spheroidal harmonic coefficients
        function obj = get_shc(obj)
            if isempty(obj.sigma)
                error("Must provide sigma in order to calculate coefficients.")
            else
                [np,nf,ns]=size(obj.sigma);
                sig = reshape(obj.sigma,np,[],1);
                sig_shc=shAna(sig);
                obj.sigma_coefficients=reshape(sig_shc,[],nf,ns);
            end
        end

        %Calculate sigma if we only have spherical/spheroidal harmonic coefficients
        function obj = get_sigma(obj)
            if isempty(obj.sigma_coefficients)
                error("Must provide sigma_coefficients in order to calculate sigma.")
            else
                [sp,nf,ns]=size(obj.sigma_coefficients);
                sig_shc = reshape(obj.sigma_coefficients,sp,[],1);
                sig=shSyn(sig_shc);
                obj.sigma=reshape(sig,[],nf,ns);
            end
        end

        %Center the (particle_number)^th spheroid at the origin and get
        %cartesian coordinates of all other targets Xt and itself, Xself.
        %If no particle number is given, simply find the coordinates for
        %all spheroids in their given positions.
        function [Xt, Xself] = get_X(obj,particle_number)
            % Get cartesian coordinates of points on spheroid surfaces, 
            % centered at particle # (particle_number).
            
            px=obj.p;
            if px == 0
                fprintf("SpheroidalParameters object was not assigned an order, p.\nUsing a default of p=8 to generate coordinates.\n");
                px=8;
            end

            ns=size(obj.centers,1); %number of particles
            np=2*px*(px+1);
            u0_arr=obj.u0;
            a_arr=obj.a;

            if length(obj.u0)==1
                u0_arr=repmat(obj.u0,ns,1);
            end
            if length(obj.a)==1
                a_arr=repmat(obj.a,ns,1);
            end
           

            if nargin == 1
                Xself=[];
                Xt=zeros(np*ns,3);
                for i=1:ns
                    if obj.oblate(i)
                        X=oblate_spheroid_shape(px,u0_arr(i),a_arr(i));
                    else
                        X=prolate_spheroid_shape(px,u0_arr(i),a_arr(i));
                    end
               
                    ci=obj.centers(i,:)';
                    if isempty(obj.Rmat)
                        thetai=obj.thetas(i);
                        phii=obj.phis(i);
    
                        Riy=[cos(thetai) 0 sin(thetai); 0 1 0; -sin(thetai) 0 cos(thetai)];
                        Riz=[cos(phii) -sin(phii) 0; sin(phii) cos(phii) 0; 0 0 1];
                        Ri=Riz*Riy;
                    else
                        Ri=obj.Rmat(:,:,i);
                    end

                    X =( Ri * X' + repmat(ci,1,np) )';
            
                    Xt((i-1)*np+1:i*np,:) = X;
                
                end

            else

                c=obj.centers(particle_number,:)';
                if isempty(obj.Rmat)
                    theta=obj.thetas(particle_number);
                    phi=obj.phis(particle_number);
                else
                    R=obj.Rmat(:,:,particle_number);
                end

                Xt=zeros((ns-1)*np,3);
                for i=1:ns
                    if obj.oblate(i)
                        X=oblate_spheroid_shape(px,u0_arr(i),a_arr(i));
                    else
                        X=prolate_spheroid_shape(px,u0_arr(i),a_arr(i));
                    end

                    if i == particle_number
                        Xself=X;
                    else
                        ci=obj.centers(i,:)';
                        if isempty(obj.Rmat)
                            thetai=obj.thetas(i);
                            phii=obj.phis(i);
    
                            Riy=[cos(thetai) 0 sin(thetai); 0 1 0; -sin(thetai) 0 cos(thetai)];
                            Riz=[cos(phii) -sin(phii) 0; sin(phii) cos(phii) 0; 0 0 1];
                            Ri=Riz*Riy;
    
                            Ry=[cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
                            Rz=[cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];
                            R=Rz*Ry;
                        else
                            Ri=obj.Rmat(:,:,i);
                        end
    
                        X = ( R' * ( Ri * X' + repmat(ci-c,1,np) ) )';
                
                        ind = i-(i>particle_number);
                        Xt((ind-1)*np+1:ind*np,:) = X;
                    end
                end
            end
        end

        %Given some input targets X, find the cartesian coordinates of
        %those targets relative to each spheroid.
        function Xt = get_X_targets(obj,X)
            ns=size(obj.centers,1); %number of particles
            np=size(X,1);
            Xt=zeros(np,3,ns);
            
            for i=1:ns

                ci=obj.centers(i,:)';
                if isempty(obj.Rmat)
                    thetai=obj.thetas(i);
                    phii=obj.phis(i);
    
                    Riy=[cos(thetai) 0 sin(thetai); 0 1 0; -sin(thetai) 0 cos(thetai)];
                    Riz=[cos(phii) -sin(phii) 0; sin(phii) cos(phii) 0; 0 0 1];
                    Ri=Riz*Riy;
                else
                    Ri=obj.Rmat(:,:,i);
                end
        
                Xt(:,:,i) = ( Ri' * (X' - repmat(ci,1,np)) )';
             
            end
        end

        %Given cell of targets X relative to each spheroid,  then find the 
        % cartesian coordinates of those targets globally.
        function Xt = set_X_targets(obj,Xcell)
            ns=size(obj.centers,1); %number of particles
            if length(Xcell)~=ns
                error('Length of input cell does not match the number of spheroids')
            end

            Xt=zeros(size(cell2mat(reshape(Xcell,[],1))));

            start=1;
            for i=1:ns
                X=Xcell{i};
                np=size(X,1);

                if np > 0

                    ci=obj.centers(i,:)';
                    if isempty(obj.Rmat)
                        thetai=obj.thetas(i);
                        phii=obj.phis(i);
    
                        Riy=[cos(thetai) 0 sin(thetai); 0 1 0; -sin(thetai) 0 cos(thetai)];
                        Riz=[cos(phii) -sin(phii) 0; sin(phii) cos(phii) 0; 0 0 1];
                        Ri=Riz*Riy;
                    else
                        Ri=obj.Rmat(:,:,i);
                    end
            
                    Xt(start:start+np-1,:) = ( Ri * X' + repmat(ci,1,np) )';
                    
                end

                start=start+np;
             
            end
        end

        %Use cutoff distance matvec_eta to split near and far spheroids.
        function separation=separate_spheroids(obj)
            % Given a cutoff distance matvec_eta, determine which particles are
            % neighbors. Store in separation matrix: -1 = self (diagonal), 0 = close,
            % 1 = far.

            if isempty(obj.matvec_eta)
                error("matvec_eta is empty. must assign a cutoff distance to categorize near and far targets.")
            elseif numel(obj.matvec_eta) >1
                error("matvec_eta should be a single value. must assign a single cutoff distance to categorize near and far targets.")
            end

            %number of points, number of spheroids
            ns=size(obj.centers,1);
            
            if 3 ~= size(obj.centers,2)
                error("centers has dimensions [~,%d] but must have 3 columns for (x,y,z) coordinates.", size(obj.centers,2))
            end
            if ns ~= size(obj.centers,1)
                error("centers is size [%d,3] but there are %d spheroids. centers must be specified for each particle.",size(obj.centers,1),ns);
            end
            % if ns~= length(obj.phis)
            %     error("phis array is length %d but there are %d spheroids. phi-rotations must be specified for each particle.",length(obj.phis),ns);
            % end
            % if ns~= length(obj.thetas)
            %     error("thetas array is length %d but there are %d spheroids. theta-rotations must be specified for each particle.",length(obj.thetas),ns);
            % end
            if ns~= size(obj.Rmat,3)
                err("Rmat array is length %d but there are %d spheroids. Rotations must be specified for each particle.",size(obj.Rmat,3),ns);
            end

            u0_arr=obj.u0;
            a_arr=obj.a;
            if length(obj.u0)==1
                u0_arr=repmat(obj.u0,1,ns);
            end
            if length(obj.a)==1
                a_arr=repmat(obj.a,1,ns);
            end

            if length(u0_arr) ~= ns
                error("u0 array is length %d but there are %d spheroids. u0 should be length 1 if all spheroids have the same scaling or length %d to specify for each particle.",length(u0_arr),ns,ns);
            end
            if length(a_arr) ~= ns
                error("'a' array is length %d but there are %d spheroids. 'a' should be length 1 if all spheroids have the same scaling or length %d to specify for each particle.",length(a_arr),ns,ns);
            end
            
            separation = -eye(ns); %-1=self, 0 =close, 1= far
            
            for i=1:ns-1
                
                %indices of all spheroids after current spheroid 'i'. Only need to
                %consider indices after 'i' since once a spheroid is considered 'self',
                %we compute all of its neighbors.
                ind=i+1:ns; 
                nt=length(ind);
                
                at = a_arr(ind); %target a's
                u0t = u0_arr(ind); %target u0's
                ct=obj.centers(ind,:); %target spheroid centers
                if isempty(obj.Rmat)
                    thetat=obj.thetas(ind); %target theta rotations
                    phit=obj.phis(ind); %target phi rotations
                end

                ai=a_arr(i); %self a
                u0i=u0_arr(i); %self a
                ci=obj.centers(i,:); %self center
                if isempty(obj.Rmat)
                    thetai=obj.thetas(i); %self theta rotation
                    phii=obj.phis(i); %self phi rotation
                    
                    % Orientation of self spheroid (parallel to axis of rotational symmetry)
                    vdi=[sin(thetai).*cos(phii) sin(thetai).*sin(phii) cos(thetai)];
                    vdi=vdi./norm(vdi);
    
                    % Orientation of target spheroids (parallel to axis of rotational symmetry)
                    vdt=[sin(thetat').*cos(phit') sin(thetat').*sin(phit') cos(thetat')];
                    vdt=vdt./repmat(vecnorm(vdt,2,2),1,3);
                else
                    Rt=obj.Rmat(:,:,ind);
                    Ri=obj.Rmat(:,:,i);
                    vdi=(Ri*[0;0;1]).';
                    vdi=vdi./norm(vdi);
                    vdt=zeros(length(ind),3);
                    for k=1:length(ind)
                        vdt(k,:)=(Rt(:,:,k)*[0;0;1]).';
                    end
                    vdt=vdt./repmat(vecnorm(vdt,2,2),1,3);
                end
            
                displacements = ct-repmat(ci,nt,1); %displacement of all targets from self
                d=vecnorm(displacements,2,2);
                cd = displacements./repmat(d,1,3);    %normalized directions of displacements
                
                % Dot product between displacement direction and orientation. 
                % If displacement is perpendicular (parallel) to orientation, then
                % approximate closest point using semi-minor (major) axis. Linear
                % interpolation to adjust for in between.
                alignment = abs(sum(vdt.*cd,2))';
                self_alignment = abs(vdi*cd');
                if ~obj.oblate
                    target_radii = at.*sqrt(u0t.^2-1) + (at.*u0t - at.*sqrt(u0t.^2-1)).*alignment;
                    self_radii = ai.*sqrt(u0i.^2-1) + (ai.*u0i - ai.*sqrt(u0i.^2-1)).*self_alignment;
                else
                    target_radii = at.*sqrt(u0t.^2+1).*alignment + (-at.*u0t + at.*sqrt(u0t.^2+1));
                    self_radii = ai.*sqrt(u0i.^2+1).*self_alignment + (-ai.*u0i + ai.*sqrt(u0i.^2+1));
                end

                closest_est = d' - self_radii - target_radii;
            
                %separate which spheroids are near (0) and which are far (1)
                sep = closest_est > obj.matvec_eta;
                separation(i,i+1:end) = sep;
                separation(i+1:end,i) = sep';
                
            end 
        end %end function separation = ...

        %Use cutoff distance matvec_eta to split near and far target points.
        function [separation,Xt]=separate_targets(obj,X)
            %{
                Given a cutoff distance matvec_eta, determine which targets are
                near. Store in separation matrix: -1 = self (diagonal), 0 = close,
                1 = far.

                The output should be
            %}

            if isempty(obj.matvec_eta)
                error("matvec_eta is empty. must assign a cutoff distance to categorize near and far targets.")
            elseif numel(obj.matvec_eta) >1
                error("matvec_eta should be a single value. must assign a single cutoff distance to categorize near and far targets.")
            end

            %number of points, number of spheroids
            ns=size(obj.centers,1);
            
            if 3 ~= size(obj.centers,2)
                error("centers has dimensions [~,%d] but must have 3 columns for (x,y,z) coordinates.", size(obj.centers,2))
            end
            if ns ~= size(obj.centers,1)
                error("centers is size [%d,3] but there are %d spheroids. centers must be specified for each particle.",size(obj.centers,1),ns);
            end
            % if ns~= length(obj.phis)
            %     error("phis array is length %d but there are %d spheroids. phi-rotations must be specified for each particle.",length(obj.phis),ns);
            % end
            % if ns~= length(obj.thetas)
            %     error("thetas array is length %d but there are %d spheroids. theta-rotations must be specified for each particle.",length(obj.thetas),ns);
            % end
            
            u0_arr=obj.u0;
            a_arr=obj.a;
            if length(obj.u0)==1
                u0_arr=repmat(obj.u0,1,ns);
            end
            if length(obj.a)==1
                a_arr=repmat(obj.a,1,ns);
            end

            if length(u0_arr) ~= ns
                error("u0 array is length %d but there are %d spheroids. u0 should be length 1 if all spheroids have the same scaling or length %d to specify for each particle.",length(u0_arr),ns,ns);
            end
            if length(a_arr) ~= ns
                error("'a' array is length %d but there are %d spheroids. 'a' should be length 1 if all spheroids have the same scaling or length %d to specify for each particle.",length(a_arr),ns,ns);
            end
            
            Xt = obj.get_X_targets(X);
            nt = size(Xt,1);
            separation = zeros(nt,ns); %0 =close, 1= far

            d=vecnorm(Xt,2,2);
            cd = Xt./repmat(d,1,3);    %normalized directions of displacements

            for i=1:ns
                ai=a_arr(i); %self a
                u0i=u0_arr(i); %self a
                
                % Dot product between displacement direction and orientation. 
                % If displacement is perpendicular (parallel) to orientation, then
                % approximate closest point using semi-minor (major) axis. Linear
                % interpolation to adjust for in between.
                alignment = cd(:,3,i);
                if ~obj.oblate
                    self_radii = ai.*sqrt(u0i.^2-1) + (ai.*u0i - ai.*sqrt(u0i.^2-1)).*alignment;
                else
                    self_radii = ai.*sqrt(u0i.^2+1).*alignment + (-ai.*u0i + ai.*sqrt(u0i.^2+1));
                end

                closest_est = d(:,:,i) - self_radii;
            
                %separate which targets are near (0) and which are far (1)
                sep = closest_est > obj.matvec_eta;
                separation(:,i) = sep;
            end 
        end

        %plot spheroids
        function plot(obj,Xt)
            X=obj.get_X;
            if nargin==1
                Xt=[];
            end
            Xplt=[X; Xt];
            colors=[repmat([0 0 1],size(X,1),1); repmat([1 0 0],size(Xt,1),1)];
            figure;
            scatter3(Xplt(:,1),Xplt(:,2),Xplt(:,3),16,colors,'filled');
            axis equal;
        end

        % Get the normal vector out of surface
        function nor = get_Norm(obj,ext_p,part_num)
            % find normal to self surface

            if nargin==1
                if obj.p == 0
                    fprintf("no p values stored, using p=8.\n ")
                    np=2*8*(8+1);
                    [theta_x,phi_x]=gl_grid(8);
                else
                    np=2*obj.p*(obj.p+1);
                    [theta_x,phi_x]=gl_grid(obj.p);
                end
                part_num = 1;
            else
                np=2*ext_p*(ext_p+1);
                [theta_x,phi_x]=gl_grid(ext_p);
            end
            ui=obj.u0(part_num);
            obl=obj.oblate(part_num);
            v_x=cos(theta_x);
            if ~obl
                nor = 1./sqrt((ui^2-1).*(ui^2-v_x.^2)).*[ui.*sqrt(ui.^2-1).*sqrt(1-v_x.^2).*cos(phi_x),...
                    ui.*sqrt(ui^2-1).*sqrt(1-v_x.^2).*sin(phi_x),...
                    (ui^2-1).*v_x];
            else
                nor = 1./sqrt((ui^2+1).*(ui^2+v_x.^2)).*[ui.*sqrt(ui.^2+1).*sqrt(1-v_x.^2).*cos(phi_x),...
                    ui.*sqrt(ui^2+1).*sqrt(1-v_x.^2).*sin(phi_x),...
                    (ui^2+1).*v_x];
            end
        end

        function [nor_rot,nor_self]=get_Norm_rot(obj,ext_p,center_particle)
            %{

            %}
            ns=size(obj.centers,1);
            np=2*ext_p*(ext_p+1);
            nor_rot=zeros(np*ns,3);
            % nor_rot = zeros(np, 3, ns);
            nor_self=[];
            if nargin<3
                % no center particle given, do self rotation only.
                for i=1:ns
                    if isempty(obj.Rmat)
                        phii=obj.phis(i); thetai=obj.thetas(i);
    
                        Riy=[cos(thetai) 0 sin(thetai); 0 1 0; -sin(thetai) 0 cos(thetai)];
                        Riz=[cos(phii) -sin(phii) 0; sin(phii) cos(phii) 0; 0 0 1];
                        Ri=Riz*Riy;
                    else
                        Ri=obj.Rmat(:,:,i);
                    end
    
                    nor_rot((i-1)*np+1:i*np,:) = (Ri * (obj.get_Norm(ext_p,i))')';
                    % nor_rot(:,:,i) = (Ri * (obj.get_Norm(ext_p,i))')';
                end
            else
                %{
                    Calculates the normal vectors of the other spheroids from the perspective
                    of the center_particle.

                    To do this, calculate the normal vectors in some global coordinate system,
                    and then apply the inverse rotation to convert into the local coordinate frame.
                %}
                for i=1:ns
                    nor_self=obj.get_Norm(ext_p,i);
                    if i==center_particle
                        nor_rot((i-1)*np+1:i*np,:)=nor_self;
                    else
                        if isempty(obj.Rmat)
                            phii=obj.phis(i); thetai=obj.thetas(i);
    
                            Riy=[cos(thetai) 0 sin(thetai); 0 1 0; -sin(thetai) 0 cos(thetai)];
                            Riz=[cos(phii) -sin(phii) 0; sin(phii) cos(phii) 0; 0 0 1];
                            Ri=Riz*Riy;
    
                            phi=obj.phis(center_particle); theta=obj.thetas(center_particle);
                            Ry=[cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
                            Rz=[cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];
                            R=Rz*Ry;
                        else
                            R=obj.Rmat(:,:,center_particle);
                            Ri=obj.Rmat(:,:,i);
                        end

                        ind = i-(i>center_particle);
                        nor_rot((ind-1)*np+1:ind*np,:) = (R' * (Ri * nor_self'))';
                        % nor_rot(:,:,ind) = (R' * (Ri * nor_self'))';
                    end
                end
            end
        end

        function nu_t = get_nu_targets(obj, nu)
            %{
                Given some input normal vectors nu, find the cartesian coordinates of
                those vectors relative to each spheroid's local frame.
            %}
            ns=size(obj.centers,1); %number of particles
            nt=size(nu,1);
            nu_t=zeros(nt,3,ns);
            
            for i=1:ns
                if isempty(obj.Rmat)
                    thetai=obj.thetas(i);
                    phii=obj.phis(i);
    
                    Riy=[cos(thetai) 0 sin(thetai); 0 1 0; -sin(thetai) 0 cos(thetai)];
                    Riz=[cos(phii) -sin(phii) 0; sin(phii) cos(phii) 0; 0 0 1];
                    Ri=Riz*Riy;
                else
                    Ri=obj.Rmat(:,:,i);
                end
        
                % Apply inverse rotation to transform from global to local frame
                nu_t(:,:,i) = (Ri' * nu')';
            end
        end

    end
end
