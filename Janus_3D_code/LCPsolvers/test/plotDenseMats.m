%%
srcFile ='amphi.lattice.n_5.p_8.cDist_2.5.allMats.mat';
res_ = load(['/Users/niru8088/scratch/Spheroidal3D-collisions/Janus_3D_code/data.11.12.2025/' srcFile]);
[~,jHi] = max(res_.ps);
[~,kHi]= min(res_.tols);
Nt = size(res_.A,1);
ttlStr = [sprintf('A($p=%d', res_.ps(jHi)) ', \epsilon_{\mathrm{gmres}}=' sprintf('%.1g)', res_.tols(kHi)) '$'];
A = res_.A;
dt = res_.dt;
[tols,kLo]= max(res_.tols);
IX = find(cellfun(@(A) ~isempty(A), res_.A(:,jHi,kHi)) & ...
cellfun(@(A) ~isempty(A), res_.A(:,5,kLo)) &...
cellfun(@(A) ~isempty(A), res_.A(:,4,kLo)) & ...
cellfun(@(A) ~isempty(A), res_.A(:,3,kLo)) & ...
cellfun(@(A) ~isempty(A), res_.A(:,2,kLo)) & ...
cellfun(@(A) ~isempty(A), res_.A(:,1,kLo))); 
I = numel(IX);
%%
ps = 3:8;
numP = numel(ps);
numTol = numel(tols);
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
    for jLo = 1:numP-1
        AHi = A{i,jHi, kHi};
        ALo = A{i,jLo, kLo};
        ALo = (ALo + ALo') /2;
        sqrtALoinv = inv(sqrtm(ALo));
        relErr(ii,jLo, kLo) = norm(AHi-ALo) / norm(AHi);
        preCond(ii,jLo, kLo) = cond(sqrtALoinv*AHi*sqrtALoinv);
        timePerHi(ii,jLo, kLo) = mean(dt{i,jHi,kHi}) / mean(dt{i,jLo,kLo});
    end
end
name = split(srcFile,'_');
name = join(name,'\_');
name = name{1};
% 
% f2 = figure;
% subpLot(1,1,1)
% heatmap(ps, tols, relErr', 'CoLorLimits', [0 1e2])
% sgtitle(sprintf('Relative Error for mc=%d', i) )
% title(name)
% set(f2, 'Position',  [0, 0, 1000, 1200])
% relErrFile = fullfile(basedir,'..','docs','fig', [srcFile(1:end-4) '_relErr.png']);
% disp(['Saving to ' relErrFile])
% saveas(f2, relErrFile);
% 
% f3 = figure;
% subpLot(1,1,1)
% heatmap(ps, tols, preCond', 'CoLorLimits', [1 1e2])
% sgtitle(['Condition  Number of $\hat{A}^{-1/2}A\hat{A}^{-1/2}$ for ' sprintf('mc=%d', i)], 'Interpreter', 'latex');
% title(name)
% set(f3, 'Position',  [0, 0, 1000, 1200])
% timePerHiFile = fullfile(basedir,'..','docs','fig', [srcFile(1:end-4) '_preCond.png']);
% disp(['Saving to ' timePerHiFile])
% saveas(f3, timePerHiFile);
% 
f4 = figure;
timePerHi = mean(timePerHi,1);
timePerHi = timePerHi(1,:,kLo);
timePerHi(end) =1;
subplot(1,1,1)
heatmap(ps, tols, timePerHi, 'CoLorLimits', [1 1e2])
sgtitle('Average Time to apply $\hat{A}$ vs $A$', 'Interpreter', 'latex');
title(name)
set(f4, 'Position',  [0, 0, 1000, 1200])
timePerHiFile = fullfile(basedir,'..','docs','fig', [srcFile(1:end-4) '_timePerHi.png']);
disp(['Saving to ' timePerHiFile])
saveas(f4, timePerHiFile);