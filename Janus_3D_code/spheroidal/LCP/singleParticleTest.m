addpath('../systems');
addpath(genpath('../../LCPsolvers'))

a = cell(1, 1);
b = cell(1,1);
c = cell(1,1);
a{1} = 1;
b{1} = 3;
c{1} = 1;
C = [0 0 3].';
r = sqrt(2)/2';
R = [1 0 0; 0 r -r ;0 r r];
rotConfig = struct('a', a, 'b', b, 'c', c,'C', C, 'R', R);

%First Try Ten Frames a second
totalSteps = 40;
timeStart = 0;
timeEnd = 4;

configurations = simulationFloor(rotConfig, totalSteps, timeStart, timeEnd);

fig = figure;
fig.Visible = 'off';
tenfps = VideoWriter('10fps.avi');
tenfps.FrameRate = 10;
tenfps.Quality = 100;
open(tenfps);

tenfpsCount = 0;
for i = 1:length(configurations)
    plotSpheroids(configurations{i}, 10);
    writeVideo(tenfps, getframe);
    [x1, x2, d] = planeEllipsoidIntersectionSigned(configurations{i}, [0 0 1].', [0 0 0].');
    if d < 0
        tenfpsCount = tenfpsCount + 1;
    end
end
close(tenfps);

configurations = simulationFloorSigned(rotConfig, totalSteps, timeStart, timeEnd);

fig = figure;
fig.Visible = 'off';
tenfpsSigned = VideoWriter('10fpsSigned.avi');
tenfpsSigned.FrameRate = 10;
tenfpsSigned.Quality = 100;
open(tenfpsSigned);

tenfpsSignedCount = 0;
for i = 1:length(configurations)
    plotSpheroids(configurations{i}, 10);
    writeVideo(tenfpsSigned, getframe);
    [x1, x2, d] = planeEllipsoidIntersectionSigned(configurations{i}, [0 0 1].', [0 0 0].');
    if d < 0
        tenfpsSignedCount = tenfpsSignedCount + 1;
    end
end
close(tenfpsSigned);

%Now try 100 frames a second

totalSteps = 400;
timeStart = 0;
timeEnd = 4;

configurations = simulationFloor(rotConfig, totalSteps, timeStart, timeEnd);

fig = figure;
fig.Visible = 'off';
hundredfps = VideoWriter('100fps.avi');
hundredfps.FrameRate = 100;
hundredfps.Quality = 100;
open(hundredfps);

hundredfpsCount = 0;
for i = 1:length(configurations)
    plotSpheroids(configurations{i}, 10);
    writeVideo(hundredfps, getframe);
    [x1, x2, d] = planeEllipsoidIntersectionSigned(configurations{i}, [0 0 1].', [0 0 0].');
    if d < 0
        hundredfpsCount = hundredfpsCount + 1;
    end
end
close(hundredfps);

configurations = simulationFloorSigned(rotConfig, totalSteps, timeStart, timeEnd);

fig = figure;
fig.Visible = 'off';
hundredfpsSigned = VideoWriter('100fpsSigned.avi');
hundredfpsSigned.FrameRate = 100;
hundredfpsSigned.Quality = 100;
open(hundredfpsSigned);

hundredfpsSignedCount = 0;
for i = 1:length(configurations)
    plotSpheroids(configurations{i}, 10);
    writeVideo(hundredfpsSigned, getframe);
    [x1, x2, d] = planeEllipsoidIntersectionSigned(configurations{i}, [0 0 1].', [0 0 0].');
    if d < 0
        hundredfpsSignedCount = hundredfpsSignedCount + 1;
    end
end
close(hundredfpsSigned);

%Now try 1000 frames a second

totalSteps = 4000;
timeStart = 0;
timeEnd = 4;

configurations = simulationFloor(rotConfig, totalSteps, timeStart, timeEnd);

fig = figure;
fig.Visible = 'off';
thousandfps = VideoWriter('1000fps.avi');
thousandfps.FrameRate = 1000;
thousandfps.Quality = 100;
open(thousandfps);

thousandfpsCount = 0;
for i = 1:length(configurations)
    plotSpheroids(configurations{i}, 10);
    writeVideo(thousandfps, getframe);
    [x1, x2, d] = planeEllipsoidIntersectionSigned(configurations{i}, [0 0 1].', [0 0 0].');
    if d < 0
        thousandfpsCount = thousandfpsCount + 1;
    end
end
close(thousandfps);

configurations = simulationFloorSigned(rotConfig, totalSteps, timeStart, timeEnd);

fig = figure;
fig.Visible = 'off';
thousandfpsSigned = VideoWriter('1000fpsSigned.avi');
thousandfpsSigned.FrameRate = 1000;
thousandfpsSigned.Quality = 100;
open(thousandfpsSigned);

thousandfpsSignedCount = 0;
for i = 1:length(configurations)
    plotSpheroids(configurations{i}, 10);
    writeVideo(thousandfpsSigned, getframe);
    [x1, x2, d] = planeEllipsoidIntersectionSigned(configurations{i}, [0 0 1].', [0 0 0].');
    if d < 0
        thousandfpsSignedCount = thousandfpsSignedCount + 1;
    end
end
close(thousandfpsSigned);


%Now try 10000 frames a second

totalSteps = 40000;
timeStart = 0;
timeEnd = 4;

configurations = simulationFloor(rotConfig, totalSteps, timeStart, timeEnd);

fig = figure;
fig.Visible = 'off';
tenthousandfps = VideoWriter('10000fps.avi');
tenthousandfps.FrameRate = 10000;
tenthousandfps.Quality = 100;
open(tenthousandfps);

tenthousandfpsCount = 0;
for i = 1:length(configurations)
    plotSpheroids(configurations{i}, 10);
    writeVideo(tenthousandfps, getframe);
    [x1, x2, d] = planeEllipsoidIntersectionSigned(configurations{i}, [0 0 1].', [0 0 0].');
    if d < 0
        tenthousandfpsCount = tenthousandfpsCount + 1;
    end
end
close(tenthousandfps);

configurations = simulationFloorSigned(rotConfig, totalSteps, timeStart, timeEnd);

fig = figure;
fig.Visible = 'off';
tenthousandfpsSigned = VideoWriter('10000fpsSigned.avi');
tenthousandfpsSigned.FrameRate = 10000;
tenthousandfpsSigned.Quality = 100;
open(thousandfpsSigned);

thousandfpsSignedCount = 0;
for i = 1:length(configurations)
    plotSpheroids(configurations{i}, 10);
    writeVideo(tenthousandfpsSigned, getframe);
    [x1, x2, d] = planeEllipsoidIntersectionSigned(configurations{i}, [0 0 1].', [0 0 0].');
    if d < 0
        tenthousandfpsSignedCount = tenthousandfpsSignedCount + 1;
    end
end
close(tenthousandfpsSigned);

%Now compute proportion of ellipsoids with negative distance to the plane 

counts = [tenfpsCount, hundredfpsCount, thousandfpsCount];
signedCounts = [tenfpsSignedCount, hundredfpsSignedCount, thousandfpsSignedCount];

totals = [40, 400, 4000];
proportions = counts./totals;
signedProportions = signedCounts./totals;

X = categorical('10', '100', '1000');
X = reordercats(X, '10', '100', '1000');
Y = [proportions.' signedProportions.'];
bar(X,Y);
saveas(gcf, 'proportions.pdf');

