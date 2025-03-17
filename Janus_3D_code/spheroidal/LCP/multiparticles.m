addpath('../systems');
addpath(genpath('../../LCPsolvers'))


initialConfiguration = RSA(10, 2, 1.5, 50);
for i =1:length(initialConfiguration)
    initialConfiguration(i).C = initialConfiguration(i).C + [0 0 5].';
end

totalSteps = 5000;
timeStart = 0;
timeEnd = 5;
configurations = simulationFloor(initialConfiguration, totalSteps, timeStart, timeEnd);


fig = figure;
fig.Visible = 'off';
hundredfpsSigned = VideoWriter('100fpsSignedMulti.avi');
hundredfpsSigned.FrameRate = 1000;
hundredfpsSigned.Quality = 100;
open(hundredfpsSigned);

hundredfpsSignedCount = 0;
for i = 1:length(configurations)
    plotSpheroids(configurations{i}, 10);
    writeVideo(hundredfpsSigned, getframe);
end
close(hundredfpsSigned);

