%single spheroid testing. The solution to the LCP has a closed form as its a constrained 1D 
% quadratic.
totalParticles = 1;
mass = 1;
%Make mobility matrix, for all the particles and the wall. The wall is modeled as having infinite mass (it doesn't move), so zeros for the corresponding parts of the mobility matrix.
mobilityMatrix = zeros(6*(totalParticles + 1), 6*(totalParticles + 1));

totalSteps = 10;
timeStart = 0;
timeEnd = 2;
timeStep = (timeEnd - timeStart)/totalSteps;


%define important values for the distance algorithms
tol = 1e-6;
maxIter = 200;
eta = 10;

addpath('../systems');
addpath(genpath('../../LCPsolvers'))
a{1} = 1;
b{1} = 3;
c{1} = 1;
C = [0 0 3].';
r = sqrt(2)/2';
R = [1 0 0; 0 r -r ;0 r r];
[R, ~] = qr(randn(3,3));
initialConfiguration = struct('a', a, 'b', b, 'c', c,'C', C, 'R', R);
configurations = cell(totalSteps,1);
velocities = cell(totalSteps,1);
quats = cell(totalSteps,1);
configurations{1} = initialConfiguration;
configuration = initialConfiguration;
velocity = zeros(6*(totalParticles + 1), 1);

quaternionConfiguration = zeros(7*(totalParticles + 1),1);
quaternionConfiguration(1:7*(totalParticles)) = rot2Quat(configuration);

for i = 1:totalParticles 
    mobilityMatrix(6*(i - 1) + 1:6*(i - 1) + 3, 6*(i - 1) + 1:6*(i - 1) + 3) = (1/mass).*(eye(3));
    mobilityMatrix(6*(i - 1) + 4:6*(i - 1) + 6, 6*(i - 1) + 4:6*(i - 1) + 6) = diag(1./(((1/5)*mass).*[initialConfiguration(i).b^2 + initialConfiguration(i).c^2 initialConfiguration(i).a^2 + initialConfiguration(i).c^2 initialConfiguration(i).a^2 + initialConfiguration(i).b^2]));
end

gravity = zeros(6*(totalParticles + 1), 1);
for i = 1:totalParticles
    gravity(6*(i - 1) + 1:6*(i - 1) + 3) = [0 0 -9.81*mass].';
end

for i = 2:totalSteps + 1
    [x1, x2, d, neighbors] = allWithWall(configuration, eta);
    D = zeros(6*2, 1);

    N = (x2{1} - x1{1})./d;
    D(1:3) = -N;

    torque = -cross(x1{1} - configuration(1).C, N);
    D(4:6) = torque;

    disp(D);
    forceExternal = zeros(6*(totalParticles + 1), 1);
    forceExternal = forceExternal + gravity;
    M = timeStep^2.*D.'*mobilityMatrix*D;
    q = d + timeStep.*D.'*velocity + timeStep^2.*D.'*mobilityMatrix*forceExternal;

    
    gamma = max([0, -q(1)/M(1)]);

    G = zeros(7*2, 6*2);
    G(1:3,1:3)  = eye(3);
    G(4:7,4:6) = constructPsi(quaternionConfiguration(4:7));

    newVelocity = velocity + timeStep.*mobilityMatrix*D*gamma + timeStep.*mobilityMatrix*forceExternal;

    newQuaternionConfiguration = quaternionConfiguration + timeStep.*G*newVelocity;
    configuration = quat2Rot(newQuaternionConfiguration(1:7*totalParticles), initialConfiguration(1).a, initialConfiguration(1).b, initialConfiguration(1).c);

    configurations{i} = configuration;
    velocities{i} = newVelocity;
    quats{i} = newQuaternionConfiguration;
    velocity = newVelocity;
    quaternionConfiguration = newQuaternionConfiguration;
    
    
end

for i = 1:length(configurations)
    clf;
    plotSpheroids(configurations{i}, 20);
    pause(0.1);
end

plotSpheroids(configurations{1}, 20);


function psi = constructPsi(q)
    P = [0 -q(4) q(3);
         q(4) 0 -q(2);
         -q(3) q(2) 0];

    psi = (1/2).*[-q(2:4).' ; 
                  q(1).*eye(3) - P];
end




