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
b = res_.b;
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
%% get warm start info
res_ = load('/Users/niru8088/scratch/Spheroidal3D-collisions/Janus_3D_code/data.11.03.2025/lcp.amphi.lattice.n_5.p_8.cDist_2.5.mat');
contactPairIX = cell(Nt,1);
for i = 1:Nt
    F = res_.lcp_list(i).F;
    if isempty(F)
        continue
    end
    n = size(F,2); % number of contact pairs
    N = size(F,1) / 6; % number of particles
    % For spheres torque is 0, so F will only be nonzero at the
    % positional places
    thesePairs = zeros(n,1);
    for ii = 1:n
        l0 = (find(F(:,ii), 1,'first')-1) / 6;
        l1 = (find(F(:,ii), 1,'last')-3) / 6;
        % linear indexing from 0
        % pair 0,1 -> 1, N,0 -> N(N-1), and so on
        thesePairs(ii) = (l0-1)*N + l1; 
    end
    contactPairIX{i} = thesePairs;
end
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
absErr = zeros(I,numP, numTol);
preCond = zeros(I,numP, numTol);
timePerHi = zeros(I,numP, numTol);
boundHolds = zeros(I,numP, numTol);
warmStartBoundHolds = zeros(I,1);
for ii = 1:I
    fprintf('- i %d\n', i)
    i = IX(ii);
    AHi = A{i,jHi, kHi};
    bHi = b{i};
    n = size(AHi,1);
    % cvx_begin quiet
    %     variable z(n) 
    %     minimize max( (zplus - zminus) .* (AHi*(zplus - zminus)) )
    %     subject to
    %         zplus+zminus == 1
    %         zplus >= 0
    %         zminus >= 0
    % cvx_end
    % cA = cvx_optval;
    if norm(AHi) > 10
         warmStartBoundHolds(ii) = NaN;
         boundHolds(ii,:,:) = NaN;
         absErr(ii,:,:) = NaN;
         preCond(ii,:,:) = NaN;
         timePerHi(ii,:,:) = NaN;
         continue
    end
    cA = min(eig(AHi));
    fprintf('- c(A) %.4g\n', cA)
    xHi = callCVX(zeros(n,1), AHi, bHi);
    if ii > 1 
        im1 = i-1; 
        A_im1 = A{im1,jHi, kHi};
        b_im1 = b{im1};
        n_im1 = numel(b_im1);
        if isempty(A_im1) || norm(A_im1) > 10 
            warmStartBoundHolds(ii) = NaN;
        else
            % Map the solution and LCP to the indicies of the smaller
            % problem
            ix_im1 = contactPairIX{i-1};
            ix_i = contactPairIX{i};
            ix_c = intersect(ix_i, ix_im1);
            n_c = numel(ix_c);
            b_ic = zeros(n_c,1);
            A_ic = zeros(n_c,n_c);
            b_im1c = zeros(n_c,1);
            A_im1c = zeros(n_c,n_c);
            for iii = 1:n_c
                jj = ix_c(iii) == ix_i;
                b_ic(iii) = bHi(jj);
                for iv = 1:n_im1
                    jv = ix_im1(iv) == ix_i;
                    A_ic(iii,iv) = AHi(jj,jv);
                end
                jj = ix_c(iii) == ix_im1;
                b_im1c(iii) = b_im1(jj);
                for iv = 1:n_im1
                    jv = ix_im1(iv) == ix_im1;
                    A_im1c(iii,iv) = A_im1(jj,jv);
                end
            end
            x_ic = callCVX(zeros(n_c,1), A_ic, b_ic);
            x_im1c = callCVX(zeros(n_c,1), A_im1c, b_im1c);
            delta = norm(A_ic-A_im1c, Inf);
            cprime = max(1, (min(eig(A_ic)) + delta)*norm(max(b_ic,0), Inf)) / (min(eig(A_ic)) - delta);
            if delta > 1e6 
                warmStartBoundHolds(ii) = NaN;
            elseif delta < min(eig(A_ic)) && ...
                    norm(x_ic - x_im1c, Inf) <= cprime * (norm(A_ic-A_im1c, Inf) + norm(b_ic-b_im1c,Inf))
                warmStartBoundHolds(ii) = 1;
            end
        end
    end
    for jLo = 1:numP
        for kLo = 1:numTol
            pLo = ps(jLo);
            tolLo = tols(kLo);
            ALo = A{i,jLo, kLo};
            delta = norm(AHi-ALo, Inf);
            if delta > 1e6 || norm(ALo) > 10
                boundHolds(ii,jLo, kLo) = NaN;
                absErr(ii,jLo, kLo) = NaN;
                preCond(ii,jLo, kLo) = NaN;
                timePerHi(ii,jLo, kLo) = NaN;
            else
                fprintf('-- p %d\n', pLo)
                fprintf('-- tol %.0e\n', tolLo)
                fprintf('-- delta %.4g\n', delta)
                xLo = callCVX(zeros(n,1), (ALo + ALo') /2, bHi);
                cprime = max(1, (cA + delta)*norm(max(bHi,0), Inf)) / (cA - delta); 
                if delta < cA && norm(xLo - xHi, Inf) <= cprime * (norm(AHi-ALo, Inf))
                    boundHolds(ii,jLo, kLo) = 1;
                end
                sqrtALoinv = inv(sqrtm(ALo));
                absErr(ii,jLo, kLo) = delta;
                kappa = cond(sqrtALoinv*AHi*sqrtALoinv);
                preCond(ii,jLo, kLo) = kappa;
                timePerHi(ii,jLo, kLo) = mean(dt{i,jHi,kHi}) / mean(dt{i,jLo,kLo});
            end
        end
    end
end
fprintf('The warm start bound holds %.2f \% of the time\n', mean(warmStartBoundHolds,1, "omitnan")*100)
boundHolds = mean(boundHolds,1, "omitnan");
boundHolds = reshape(boundHolds, numP, numTol);
absErr = mean(absErr,1, "omitnan");
absErr = reshape(absErr, numP, numTol);
preCond = mean(preCond,1, "omitnan");
preCond = reshape(preCond, numP, numTol);
timePerHi = mean(timePerHi,1, "omitnan");
timePerHi = reshape(timePerHi, numP, numTol);
%%
name = split(srcFile,'_');
name = join(name,'\_');
name = name{1};
%% 
f1 = figure;
subplot(1,1,1)
h = heatmap(ps, tols, 100*boundHolds');
h.CellLabelFormat = '%.0f%%';
colormap(viridis)
clim([0 100]);
sgtitle('$\|\mathbf{x} - \hat{\mathbf{x}}\|_\infty \leq c^\prime \|\mathbf{A} - \hat{\mathbf{A}}\|_\infty$', 'Interpreter', 'latex')
% title(name)
set(f1, 'Position',  [0, 0, 1000, 1200])
xlabel('p')
ylabel('\epsilon')
% set(gca,'Interpreter','latex')
fontsize(f1, 30, 'points')
absErrFile = fullfile(basedir,'..','docs','fig', [srcFile(1:end-4) '_boundsHold.png']);
disp(['Saving to ' absErrFile])
saveas(f1, absErrFile);
%% 
f2 = figure;
subplot(1,1,1)
h = heatmap(ps, tols, 100*absErr');
h.CellLabelFormat = '%.2g%%';
colormap(viridis)
clim([0 5]);
sgtitle('Absolute Error $\|\mathbf{A} - \hat{\mathbf{A}}\|_\infty$', 'Interpreter', 'latex')
% title(name)
set(f2, 'Position',  [0, 0, 1000, 1200])
xlabel('p')
ylabel('\epsilon')
% set(gca,'Interpreter','latex')
fontsize(f2, 30, 'points')
absErrFile = fullfile(basedir,'..','docs','fig', [srcFile(1:end-4) '_absErr.png']);
disp(['Saving to ' absErrFile])
saveas(f2, absErrFile);
%%
f3 = figure;
subplot(1,1,1)
h = heatmap(ps, tols, preCond');
h.CellLabelFormat = '%.2g';
colormap(viridis)
clim([1 2.5]);
sgtitle('Condition  Number of $\hat{A}^{-1/2}A\hat{A}^{-1/2}$', 'Interpreter', 'latex');
% title(name)
set(f3, 'Position',  [0, 0, 1000, 1200])
xlabel('p')
ylabel('\epsilon')
% set(gca,'Interpreter','latex')
fontsize(f3, 30, 'points')
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
sgtitle('Average Time to apply $\hat{\mathbf{A}}$ vs $\mathbf{A}$', 'Interpreter', 'latex');
% title(name)
set(f4, 'Position',  [0, 0, 1000, 1200])
xlabel('p')
ylabel('\epsilon')
% set(gca,'Interpreter','latex')
fontsize(f4, 30, 'points')
timePerHiFile = fullfile(basedir,'..','docs','fig', [srcFile(1:end-4) '_timePerHi.png']);
disp(['Saving to ' timePerHiFile])
saveas(f4, timePerHiFile);