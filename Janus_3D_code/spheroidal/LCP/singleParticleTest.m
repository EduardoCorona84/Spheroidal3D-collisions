a = cell(1, 1);
b = cell(1,1);
c = cell(1,1);
a{1} = 1;
b{1} = 1;
c{1} = 2;

addpath('../systems');
addpath(genpath('../../LCPsolvers'))

a{1} = 1;
b{1} = 3;
c{1} = 1;
C = [0 0 3].';
r = sqrt(2)/2';
R = [1 0 0; 0 r -r ;0 r r];
initialConfiguration = struct('a', a, 'b', b, 'c', c,'C', C, 'R', R);


%C = [0 0 5].';
%[R, ~] = qr(randn(3,3));
rotConfig = struct('a', a, 'b', b, 'c', c,'C', C, 'R', R);

totalSteps = 100;
timeStart = 0;
timeEnd = 2;

configurations = simulation(rotConfig, totalSteps, timeStart, timeEnd);

for i = 1:length(configurations)
    clf;
    plotSpheroids(configurations{i}, 20);
    pause(0.1);
end