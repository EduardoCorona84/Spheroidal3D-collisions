function configurations = simulationFloorSigned(initialConfiguration, totalSteps, timeStart, timeEnd, initialVelocity, density, dragType)

    if nargin < 7
        dragType = 'none';
    end

    if nargin < 6
        density = 1.0;
    end

    %add paths to needed functions
    addpath('../systems');
    addpath(genpath('../../LCPsolvers'))

    %define important values for the problem
    timeStep = (timeEnd - timeStart)/totalSteps;
    totalParticles = length(initialConfiguration);
    planeNormal = [0 0 1].';
    planePoint = [0 0 0].';

    %initialize variables
    configurations = cell(totalSteps,1);
    mobilityMatrix = zeros(6*(totalParticles + 1), 6*(totalParticles + 1));
    gravity = zeros(6*(totalParticles + 1), 1);

    %fill in first initial values
    configurations{1} = initialConfiguration;
    configuration = initialConfiguration;

    if nargin < 5
        velocity = zeros(6*(totalParticles + 1), 1);
    else
        velocity = initialVelocity;
    end

    quaternionConfiguration = zeros(7*(totalParticles + 1),1);
    quaternionConfiguration(1:7*(totalParticles)) = rot2Quat(configuration);

    %define important values for the distance algorithms
    tol = 1e-13;
    maxIter = 2000;
    eta = 1;

    %fill in the mobility matrix and gravity force vector

    for i = 1:totalParticles 
        mass = density*4/3*pi*initialConfiguration(i).a*initialConfiguration(i).b*initialConfiguration(i).c;
        mobilityMatrix(6*(i - 1) + 1:6*(i - 1) + 3, 6*(i - 1) + 1:6*(i - 1) + 3) = (1/mass).*(eye(3));
        mobilityMatrix(6*(i - 1) + 4:6*(i - 1) + 6, 6*(i - 1) + 4:6*(i - 1) + 6) = diag(((mass/5).*[initialConfiguration(i).b^2 + initialConfiguration(i).c^2 initialConfiguration(i).a^2 + initialConfiguration(i).c^2 initialConfiguration(i).a^2 + initialConfiguration(i).b^2]).^(-1));

        gravity(6*(i - 1) + 1:6*(i - 1) + 3) = [0 0 -9.81*mass].';
    end

    %iterate through the time steps.
    %i = 1 corresponds to the initial condition (t = 0), i = 2 corresponds to the t = 0 + timeStep, and so on.
    for i = 2:totalSteps + 1

        %find needed info for pairs of particles
        [x1Pairs, x2Pairs, distancesPairs, neighborsPairs] = allWithPairs(configuration, tol, maxIter, eta, 'mb');
        totalPairs = nnz(neighborsPairs)./2;
        
        %find needed info for particles and the wall
        
        [x1Wall, x2Wall, distancesWall, neighborsWall] = allWithWallSigned(configuration, planeNormal, planePoint, eta);
        totalWall = nnz(neighborsWall);
        
        %compile the previous information into combined arrays for both pairs and the wall. The wall is treated as if its the last particle (in the ordering).
        x1 = cell(totalParticles + 1, totalParticles + 1);
        x2 = cell(totalParticles + 1, totalParticles + 1);
        distances = zeros(totalParticles + 1, totalParticles + 1);
        neighbors = logical(false(totalParticles + 1, totalParticles + 1));

        x1(1:totalParticles, 1:totalParticles) = x1Pairs;
        x1(1:totalParticles, totalParticles + 1) = x1Wall;
        x1(totalParticles + 1, 1:totalParticles) = x2Wall.';

        x2(1:totalParticles, 1:totalParticles) = x2Pairs;
        x2(1:totalParticles, totalParticles + 1) = x2Wall;
        x2(totalParticles + 1, 1:totalParticles) = x1Wall.';

        distances(1:totalParticles, 1:totalParticles) = distancesPairs;
        distances(1:totalParticles, totalParticles + 1) = distancesWall;
        distances(totalParticles + 1, 1:totalParticles) = distancesWall.';

        neighbors(1:totalParticles, 1:totalParticles) = neighborsPairs;
        neighbors(1:totalParticles, totalParticles + 1) = neighborsWall;
        neighbors(totalParticles + 1, 1:totalParticles) = neighborsWall.';

        %Initialize the D matrix, this matrix maps impulses to forces and torques. We only consider pairs that are close to contact. 
        D = zeros(6*(totalParticles + 1), totalPairs + totalWall);

        %We now want to fill the D matrix. The ordering of the pairs is as follows:
        %Start with the 1st particle, then consider all the pairs involving the first particle of the form (1,x) where x from smallest to largest x. Repeat this process with the 2nd particle, this time starting x from greater than or equal to 3 (to not double count the (1,2)) pair
        currentPair = 1;
        contactDistances = zeros(totalPairs + totalWall,1);
        for j = 1:totalParticles

            %find all the neighbors of particle j
            localNeighbors = find(neighbors(j,j + 1:totalParticles + 1)) + j ;

            for k = 1:length(localNeighbors)
                if localNeighbors(k) == totalParticles + 1
                    N = -[0 0 1].';

                    %find the local coordinates, ignoring the wall
                    localCoordinates1 = x1{j,localNeighbors(k)} - configuration(j).C;

                    %Find the torque on each body, in the body rotational frame, ignoring the wall
                    torque1 = configuration(j).R.'*(-cross(localCoordinates1, N));
    
                    %Fill in Jacobian Matrix values
                    %Force and torque on body 1
                    D(6*(j - 1) + 1:6*(j - 1) + 3, currentPair) = -N;
                    D(6*(j - 1) + 4:6*(j - 1) + 6, currentPair) = torque1;

                    contactDistances(currentPair) = distances(j,localNeighbors(k));
                    currentPair = currentPair + 1;

                else
                    if distances(j,localNeighbors(k)) < tol
                        N = computeNormal(configuration(j), x2{j, localNeighbors(k)});
                    else
                        N = (x2{j,localNeighbors(k)} - x1{j,localNeighbors(k)})./distances(j,localNeighbors(k));
                    end
  
                    %Find local coordinates of contact points in the global rotational frame.
                    localCoordinates1 = x1{j,localNeighbors(k)} - configuration(j).C;
                    localCoordinates2 = x2{j,localNeighbors(k)} - configuration(localNeighbors(k)).C;
        
    
                    %Find the torque on each body, in the body rotational frame.
                    torque1 = configuration(j).R.'*(-cross(localCoordinates1, N));
                    torque2 = configuration(localNeighbors(k)).R.'*(cross(localCoordinates2, N));
    
                    %Fill in Jacobian Matrix values
                    %Force and torque on body 1
                    D(6*(j - 1) + 1:6*(j - 1) + 3, currentPair) = -N;
                    D(6*(j - 1) + 4:6*(j - 1) + 6, currentPair) = torque1;
    
                    %Force and torque on body 2
                    D(6*(localNeighbors(k) - 1) + 1:6*(localNeighbors(k) - 1) + 3, currentPair) = N;
                    D(6*(localNeighbors(k) - 1) + 4:6*(localNeighbors(k) - 1) + 6, currentPair) = torque2;
    
                    %Iterate to the next pair
                    contactDistances(currentPair) = distances(j,localNeighbors(k));
                    currentPair = currentPair + 1;
                end
            end

        end


        %Define and fill the external force, starting with the torques seen as "fictious" forces as we are in the rotating local frame.
        forceExternal = zeros(6*(totalParticles + 1), 1);

        for j = 1:totalParticles
            forceExternal(6*(j - 1) + 4:6*(j - 1) + 6) = -cross(velocity(6*(j - 1) + 4:6*(j - 1) + 6),diag(mobilityMatrix(6*(j - 1) + 4:6*(j - 1) + 6,6*(j - 1) + 4:6*(j - 1) + 6)).^(-1).*velocity(6*(j - 1) + 4:6*(j - 1) + 6));
        end

        %add drag forces to the external forces
        for j = 1:totalParticles
            forceExternal(6*(j - 1) + 1:6*(j - 1) + 3) = computeDrag(configuration(j), velocity(6*(j - 1) + 1:6*(j - 1) + 3), dragType);
        end
        
        %add gravity
        forceExternal = forceExternal + gravity;
        M = timeStep^2.*D.'*mobilityMatrix*D;
        q = contactDistances + timeStep.*D.'*velocity + timeStep^2.*D.'*mobilityMatrix*forceExternal;

        %Now solve the LCP 
        [gamma, ~, ~, ~, ~, ~] = BBPGD(M, q, zeros(size(q)), 2000, 1e-8);

        %Construct the G matrix
        G = zeros(7*(totalParticles + 1), 6*(totalParticles + 1));
        for j = 1:totalParticles
            G(7*(j - 1) + 1:7*(j - 1) + 3, 6*(j - 1) + 1:6*(j - 1) + 3) = eye(3);
            G(7*(j - 1) + 4:7*(j - 1) + 7, 6*(j - 1) + 4:6*(j - 1) + 6) = constructPsiLocal(quaternionConfiguration(7*(j - 1) + 4:7*(j - 1) + 7));
        end

        %update quantities

        newVelocity = velocity + timeStep.*mobilityMatrix*D*gamma + timeStep.*mobilityMatrix*forceExternal;
        newQuaternionConfiguration = quaternionConfiguration + timeStep.*G*newVelocity;
        configuration = quat2Rot(newQuaternionConfiguration(1:7*totalParticles), initialConfiguration(1).a, initialConfiguration(1).b, initialConfiguration(1).c); 

        configurations{i} = configuration;
        velocity = newVelocity;
        quaternionConfiguration = newQuaternionConfiguration;
    end


end
