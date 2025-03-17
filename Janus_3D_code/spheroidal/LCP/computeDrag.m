function drag = computeDrag(par, velocity, type)

    %Find cross sectional areas
    area = computeCrossSectionalArea(par);

    if strcmp(type,'linear')
        %Linear drag, assuming in water with dynamic viscosity of 0.0010518 Pa*s (of water at 18 degrees Celsius)
        viscosity = 0.0010518;
        Deq = 2*(par.a*par.b*par.c)^(1/3);
        drag = -(viscosity/(2*Deq)).*area.*velocity;
    elseif strcmp(type,'quadratic')
        %Quadratic drag
        %Ideally, we would have some sort of shape parameter that depended on the orientation of the particle. I have not found a good way to do this yet, so we will leave that for now.
        %assuming density of air
        directions = sign(velocity);
        density = 1.225;
        drag = -0.5*density*area.*velocity.^2.*directions;
    else 
        drag = zeros(3,1);
    end

end

function area = computeCrossSectionalArea(par)
    %This function computes the cross sectional area of the ellipsoid projected onto each of the 3 coordinate planes. 
    area = zeros(3,1);
    A = par.R*diag([par.a par.b par.c].^(-2))*par.R.';
    for i = 1:3
        area(i) = sqrt(A(i,i))/sqrt(det(A));
    end
    area = pi.*area;
end