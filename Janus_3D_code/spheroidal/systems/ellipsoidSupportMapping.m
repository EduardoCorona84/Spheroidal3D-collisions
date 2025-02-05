function w = ellipsoidSupportMapping(v, a, b, c, C, R)
    %v is the vector getting mapped via the support mapping.
    %a, b, c, are the semi axes of the ellipsoid
    %C and R are the center and rotation matrix 

    vRotated = R'*v;
    D = diag([a b c]);
    wTranslated = ((D.^2)*vRotated)/(norm(D*vRotated));
    w = R*wTranslated + C;
    
end