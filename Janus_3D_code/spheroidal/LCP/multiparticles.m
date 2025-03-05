addpath('../systems');
addpath(genpath('../../LCPsolvers'))
spheroids = RSA(10, 2, 2, 4);

for i = 1:length(spheroids)
    spheroids(i).C = spheroids(i).C + [0 0 5].';
end


totalSteps  = 1000;
timeStart = 0;
timeEnd = 2;

configurations = simulation(spheroids, totalSteps, timeStart, timeEnd);

for i = 1:length(configurations)
    clf;
    plotSpheroids(configurations{i}, 15);
    pause(0.001);
end