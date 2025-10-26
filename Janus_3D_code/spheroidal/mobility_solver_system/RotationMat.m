function M = RotationMat(wh,t)
    %{
    An implementation of Rodrigues' rotation matrix formula.

    Inputs
    wh - (double) angular velocity vector
    t  - (double) timestep

    Outputs
    M - (double) 3x3 rotation matrix
    %}
    nwh = norm(wh); 
    t = nwh*t; 
    wh = wh./nwh; 
    
    M = [
        1-(wh(2)^2+wh(3)^2)*(1-cos(t)) , wh(2)*wh(1)*(1-cos(t))-wh(3)*sin(t) , wh(1)*wh(3)*(1-cos(t))+wh(2)*sin(t);...
        wh(1)*wh(2)*(1-cos(t))+wh(3)*sin(t),1-(wh(1)^2+wh(3)^2)*(1-cos(t)),wh(2)*wh(3)*(1-cos(t))-wh(1)*sin(t);...
        wh(1)*wh(3)*(1-cos(t))-wh(2)*sin(t),wh(2)*wh(3)*(1-cos(t))+wh(1)*sin(t),1-(wh(2)^2+wh(1)^2)*(1-cos(t))
    ];
end