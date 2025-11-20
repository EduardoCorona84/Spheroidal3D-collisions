%%
[dirname, basedir] = setPaths();
srcFile ='amphi.lattice.n_5.p_8.cDist_2.5.allMats.mat';
res_ = load(['/Users/niru8088/scratch/Spheroidal3D-collisions/Janus_3D_code/goodData/' srcFile]);
[~,jHi] = max(res_.ps);
[~,kHi]= min(res_.tols);
Nt = size(res_.A,1);
ttlStr = [sprintf('A($p=%d', res_.ps(jHi)) ', \epsilon_{\mathrm{gmres}}=' sprintf('%.1g)', res_.tols(kHi)) '$'];
ps = res_.ps;
tols = res_.tols;
A = res_.A;
dt = res_.dt;
numP = numel(ps);
numTol = numel(tols);
% find the time steps where we have all dense A mats
mask = true(Nt,1); 
for j = 1:numP 
    for k = 1:numTol
        mask = mask & cellfun(@(A) ~isempty(A), res_.A(:,j,k)) ; 
    end
end
IX = find(mask);
I = numel(IX);
%% Timings pLot
% figure()
% f1 = figure;
% for j = 1:numP
%     for k = 1:numTol
%         subpLot(numP, numTol, numTol*(j-1) + k);
%         edges = 10.^(-6:.5:3);
%         [counts,edges] = Histcounts(dt{i,j,k},edges);
%         g = Histogram('BinEdges',edges,'BinCounts',counts);
%         set(gca, "Xscale", "Log")
%         xticks(edges)
%         title({ ...
%             sprintf('p = %d, tol = %.0e', ps(j), tols(k)), ...
%             sprintf('mean = %.3f',mean(dt{i,j,k})),...
%             sprintf('std = %.3f', std(dt{i,j,k}))}...
%         );
%     end
% end
% name = split(srcFile,'_');
% name = join(name,'\_');
% name = name{1};
% sgtitle({sprintf('Run Time for mc=%d', i),name})
% set(f1, 'Position',  [0, 0, 1000, 1200])
% timingFile = fullfile(basedir,'..','docs','fig', [srcFile(1:end-4) '_runTime.png']);
% disp(['Saving to ' timingFile])
% saveas(f1, timingFile);
% Relative Error
relErr = zeros(I,numP, numTol);
preCond = zeros(I,numP, numTol);
timePerHi = zeros(I,numP, numTol);
for ii = 1:I
    i = IX(ii);
    for jLo = 1:numP
        for kLo = 1:numTol
            pLo = ps(jLo);
            tolLo = tols(kLo);
            AHi = A{i,jHi, kHi};
            ALo = A{i,jLo, kLo};
            % ALo = (ALo + ALo') /2;
            re = norm(AHi-ALo) / norm(AHi);
            if re > 1e6 || norm(ALo) > 10 || norm(AHi) > 10
                relErr(ii,jLo, kLo) = NaN;
                preCond(ii,jLo, kLo) = NaN;
                timePerHi(ii,jLo, kLo) = NaN;
            else
                sqrtALoinv = inv(sqrtm(ALo));
                relErr(ii,jLo, kLo) = re;
                kappa = cond(sqrtALoinv*AHi*sqrtALoinv);
                preCond(ii,jLo, kLo) = kappa;
                timePerHi(ii,jLo, kLo) = mean(dt{i,jHi,kHi}) / mean(dt{i,jLo,kLo});
            end
        end
    end
end
relErr = mean(relErr,1, "omitnan");
relErr = reshape(relErr, numP, numTol);
preCond = mean(preCond,1, "omitnan");
preCond = reshape(preCond, numP, numTol);
timePerHi = mean(timePerHi,1, "omitnan");
timePerHi = reshape(timePerHi, numP, numTol);
%%
name = split(srcFile,'_');
name = join(name,'\_');
name = name{1};
%% 
f2 = figure;
subplot(1,1,1)
h = heatmap(ps, tols, 100*relErr');
h.CellLabelFormat = '%.0f%%';
colormap(viridis)
clim([0 100]);
sgtitle('Relative Error')
title(name)
set(f2, 'Position',  [0, 0, 1000, 1200])
relErrFile = fullfile(basedir,'..','docs','fig', [srcFile(1:end-4) '_relErr.png']);
disp(['Saving to ' relErrFile])
saveas(f2, relErrFile);
%%
f3 = figure;
subplot(1,1,1)
h = heatmap(ps, tols, preCond');
h.CellLabelFormat = '%.2f';
colormap(viridis)
clim([1 2.5]);
sgtitle('Condition  Number of $\hat{A}^{-1/2}A\hat{A}^{-1/2}$', 'Interpreter', 'latex');
title(name)
set(f3, 'Position',  [0, 0, 1000, 1200])
timePerHiFile = fullfile(basedir,'..','docs','fig', [srcFile(1:end-4) '_preCond.png']);
disp(['Saving to ' timePerHiFile])
saveas(f3, timePerHiFile);
%%
f4 = figure;
subplot(1,1,1)
h=heatmap(ps, tols, timePerHi');
h.CellLabelFormat = '%.2f';
colormap(viridis)
clim([1 50])
sgtitle('Average Time to apply $\hat{A}$ vs $A$', 'Interpreter', 'latex');
title(name)
set(f4, 'Position',  [0, 0, 1000, 1200])
timePerHiFile = fullfile(basedir,'..','docs','fig', [srcFile(1:end-4) '_timePerHi.png']);
disp(['Saving to ' timePerHiFile])
saveas(f4, timePerHiFile);