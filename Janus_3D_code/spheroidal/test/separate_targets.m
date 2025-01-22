function separation=separate_targets(obj,matvec_eta)
% Given a cutoff distance matvec_eta, determine which particles are
% neighbors. Store in separation matrix: -1 = self (diagonal), 0 = close,
% 1 = far.

%number of points, number of spheroids
[np,ns]=size(obj.sigma);

u0_arr=obj.u0;
a_arr=obj.a;
if length(obj.u0)==1
    u0_arr=repmat(obj.u0,ns,1);
end
if length(obj.a)==1
    a_arr=repmat(obj.a,ns,1);
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
    ct=obj.centers(ind,:); %target centers
    % thetat=obj.thetas(ind); %target theta rotations
    % phit=obj.phis(ind); %target phi rotations
    
    ai=obj.a(i); %self a
    u0i=obj.u0(i); %self a
    ci=obj.centers(i,:); %self center
    % thetai=obj.thetas(i); %self theta rotation
    % phii=obj.phis(i); %self phi rotation
    
    % % Orientation of self spheroid (parallel to axis of rotational symmetry)
    % vdi=[sin(thetai).*cos(phii) sin(thetai).*sin(phii) cos(thetai)];
    % vdi=vdi./norm(vdi);
    % 
    % % Orientation of target spheroids (parallel to axis of rotational symmetry)
    % vdt=[sin(thetat').*cos(phit') sin(thetat').*sin(phit') cos(thetat')];
    % vdt=vdt./repmat(vecnorm(vdt,2,2),1,3);

    Rt=obj.Rmat(:,:,ind);
    Ri=obj.Rmat(:,:,i);
    vdi=(Ri*[0;0;1]).';
    vdi=vdi./norm(vdi);
    vdt=zeros(length(ind),3);
    for k=1:length(ind)
        vdt(k,:)=(Rt(:,:,k)*[0;0;1]).';
    end
    vdt=vdt./repmat(vecnorm(vdt,2,2),1,3);

    displacements = ct-repmat(ci,nt,1); %displacement of all targets from self
    d=vecnorm(displacements,2,2);
    cd = displacements./repmat(d,1,3);    %normalized directions of displacements
    
    % Dot product between displacement direction and orientation. 
    % If displacement is perpendicular (parallel) to orientation, then
    % approximate closest point using semi-minor (major) axis. Linear
    % interpolation to adjust for in between.
    alignment = abs(sum(vdt.*cd,2))';
    self_alignment = abs(vdi*cd');
    target_radii = at.*sqrt(u0t.^2-1) + (at.*u0t - at.*sqrt(u0t.^2-1)).*alignment;
    self_radii = ai.*sqrt(u0i.^2-1) + (ai.*u0i - ai.*sqrt(u0i.^2-1)).*self_alignment;
    
    closest_est = d' - self_radii - target_radii;

    %separate which targets are near (0) and which are far (1)
    sep = closest_est > matvec_eta;
    separation(i,i+1:end) = sep;
    separation(i+1:end,i) = sep';
    
end

end