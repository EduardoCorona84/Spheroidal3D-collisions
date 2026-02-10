function [C,B,D,A,L] = Build_SpheroidalAuxMats(Wg,Xg,Xc,np,n3)
%{
Builds the rigid-body motion completion operator. Note that these operators
assume the input density is in an interleaved format; i.e. for a body with
density sigma, we have

[sigma_x1, sigma_y1, sigma_z1, sigma_x2, ...]

The operator of interest is L, which projects a surface density
field onto the space of rigid-body motions. For the i-th particle, this
operator is defined as,

L_i[mu_i](x) = (1/|Gamma_i|) * \int(mu_i(y) dS_y) ...
             + tau_i^-1 * (\int((y-x_i^c) x mu_i(y) dS_y)) x (x-x_i^c)

Note that the moment-of-inertia tensor should always be diagonal since in
the BIE, we assume that the source-to-source block is in the local frame of
that source spheroid (so the off-diagonal terms are zero in the tensor).

Inputs:
Wg : np*n3 x 1 array
    quadrature weights associated with each spheroid
Xg : np*n3 x 3 array
    discretization points
Xc : n3 x 3 array
    center(s) of spheres
np : scalar
    number of discretization points per spheroid
n3 : scalar
    number of spheroids

Outputs:
C : 6*n3 x 3*np*n3 array
    maps a surface density vector 'sigma' to the total force and torque on 
    each spheroid; i.e. C*sigma = [F_1; T_1; F_2; T_2; ...], where F_i is
    the total force and T_i is the total torque on spheroid i.
B : 6*n3 x 3*np*n3 array
    the operator of interest is the transpose B^T, which maps a vector of
    forces and torques to a rigid-body velocity field on the surfaces of
    the spheroids.
D : 6*n3 x 3*np*n3 array
    unweighted version of C
A : 6*n3 x 3*np*n3 array
    weighted version of B
L : 3*np*n3 x 3*np*n3 array
    the projection operator L = B^T * C. projects any surface density onto
    the space of densities that correspond to rigid body motion
%}

N = 3*np*n3; 
% I = row index, J = column index, V = value
IC = []; JC = []; VC = []; 
VD = []; VB = []; VA = [];

if nargout==5
    IL = []; JL = []; VL = []; 
end

for i=1:n3
    % Select indices for current body.
    idx = (1:np)+np*(i-1); 
    X = Xg(idx,:);
    W = Wg(idx); 

    % Get the surface area of current spheroid; i.e. |Gamma_i|
    oW = ones(size(W));
    sumW = sum(W); 
        
    if ~isempty(Xc)
        % Center the current body at the origin (i.e. place it in the local frame)
        X = X - repmat(Xc(i,:),np,1); 
    end

    % Get the global column indices for the x,y,z components of the density
    % vector for the i-th spheroid.
    % Note the ordering:
    % [sig_x1, sig_y1, sig_z1, sig_x2, ...]
    indx = (1:3:3*np)+3*np*(i-1); 
    indy = (2:3:3*np)+3*np*(i-1); 
    indz = (3:3:3*np)+3*np*(i-1);

    % C*sigma = [int{sigma} ; int{X \times sigma}]   
    ICloc = [reshape(repmat(6*(i-1)+(1:6),np,1),[],1) ; reshape(repmat(6*(i-1)+(4:6),np,1),[],1)];
    IC = [IC ; ICloc];
    JCloc = [indx'; indy'; indz'; ...
             indy'; indz'; indx'; ...
             indz'; indx'; indy'
    ];
    JC = [JC ; JCloc];

    VCloc = [W; W; W; ...
          -W.*X(:,3); -W.*X(:,1); -W.*X(:,2); ...
          W.*X(:,2);  W.*X(:,3); W.*X(:,1)
    ];
    VC = [VC ; VCloc];

    % First three vectors are just integrals of f for each coordinate
    VD = [VD;
          oW; oW; oW; ...
          -X(:,3); -X(:,1); -X(:,2); ...
          X(:,2);  X(:,3); X(:,1)
    ];

    tau_xx = sum(W.*X(:,3).^2)+sum(W.*X(:,2).^2); 
    tau_yy = sum(W.*X(:,1).^2)+sum(W.*X(:,3).^2); 
    tau_zz = sum(W.*X(:,1).^2)+sum(W.*X(:,2).^2); 

    
    VBloc = [oW./sumW; oW./sumW; oW./sumW ; ...
          -X(:,3)./tau_xx; -X(:,1)./tau_yy; -X(:,2)./tau_zz; ...
          X(:,2)./tau_xx;  X(:,3)./tau_yy; X(:,1)./tau_zz
    ];
    VB = [VB ; VBloc];

    Wav = (1/sum(W))*W; 
    VA = [VA; Wav; Wav; Wav ; ...
         -(1/tau_xx)*W.*X(:,3); -(1/tau_yy)*W.*X(:,1); -(1/tau_zz)*W.*X(:,2);...
         (1/tau_xx)*W.*X(:,2);  (1/tau_yy)*W.*X(:,3); (1/tau_zz)*W.*X(:,1)];

    % We divide by the integral of ((X-X_c) (x) 1)^2 over Sc. 
    idxM = (1:3*np)+3*np*(i-1); 
    Cloc = zeros(6,3*np); Bloc = zeros(6,3*np); 

    Cloc((ICloc-6*(i-1)) + 6*(JCloc-3*np*(i-1)-1)) = VCloc; 
    Bloc((ICloc-6*(i-1)) + 6*(JCloc-3*np*(i-1)-1)) = VBloc;

    if nargout==5 
        [JJ,II] = meshgrid(idxM,idxM); 
        IL = [IL ; II(:)]; 
        JL = [JL ; JJ(:)]; 
        Lb = Bloc'*Cloc; % B^T*C
        VL = [VL ; Lb(:)];
    end
end

C = sparse(IC,JC,VC,6*n3,N); 
B = sparse(IC,JC,VB,6*n3,N);
D = sparse(IC,JC,VD,6*n3,N);
A = sparse(IC,JC,VA,6*n3,N); 

if nargout==5
    L = sparse(IL,JL,VL,N,N);
end

end