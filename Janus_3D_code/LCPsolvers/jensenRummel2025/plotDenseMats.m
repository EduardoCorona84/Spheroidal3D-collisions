function plotDenseMats(root, prefix)
[dirname, basedir] = setPaths();
srcFile = fullfile(root, [prefix '.denseMats.allMats.mat']);
res_ = load(srcFile);
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
res_ = load(fullfile(root, [prefix '.mat']));
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

absErr = zeros(I,numP, numTol);
preCond = zeros(I,numP, numTol);
timePerHi = zeros(I,numP, numTol);
boundHolds = zeros(I,numP, numTol);
warmStartBoundHolds = zeros(I,1);
lcpOpts.solver = 'proxquasinewton';
lcpOpts = defaultLCPOpts(lcpOpts);
for ii = 1:I
    fprintf('- i %d\n', i)
    i = IX(ii);
    AHi = A{i,jHi, kHi};
    bHi = b{i};
    n = size(AHi,1);
    tic
    cA = getC_A(AHi);
    fprintf('- cA computation time: %.4g sec\n', toc);
    if norm(AHi) > 10 || isempty(AHi)
        warmStartBoundHolds(ii) = NaN;
        boundHolds(ii,:,:) = NaN;
        absErr(ii,:,:) = NaN;
        preCond(ii,:,:) = NaN;
        timePerHi(ii,:,:) = NaN;
        continue
    end
    %% We know we have this eigen value bound
    fprintf('-- lambda_min/n <= cA <= lambda_min\n');
    fprintf('-- %.4g <= %.4g <= %.4g\n',  min(eig(AHi)) /n , cA ,min(eig(AHi)));
    tic
    lcpOpts.b = bHi; lcpOpts.A = @(x) AHi*x;
    fg = @(x, Ax) quadraticLoss(x, AHi, bHi, Ax);
    xHi = proxQuasiNewton(fg, zeros(n,1),lcpOpts);
    fprintf('- xHi computation time: %.4g sec\n', toc);
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
                for iv = 1:n_c
                    jv = ix_c(iv) == ix_i;
                    A_ic(iii,iv) = AHi(jj,jv);
                end
                jj = ix_c(iii) == ix_im1;
                b_im1c(iii) = b_im1(jj);
                for iv = 1:n_c
                    jv = ix_c(iv) == ix_im1;
                    A_im1c(iii,iv) = A_im1(jj,jv);
                end
            end
            % x_ic = callCVX(zeros(n_c,1), A_ic, b_ic);
            % x_im1c = callCVX(zeros(n_c,1), A_im1c, b_im1c);
            lcpOpts.b = b_ic; lcpOpts.A = @(x) A_ic*x;
            fg_ic = @(x, Ax) quadraticLoss(x, A_ic, b_ic, Ax);
            x_ic = proxQuasiNewton(fg_ic, zeros(n_c,1), lcpOpts);
            lcpOpts.b = b_im1c; lcpOpts.A = @(x) A_im1c*x;
            fg_im1c = @(x, Ax) quadraticLoss(x, A_im1c, b_im1c, Ax);
            x_im1c = proxQuasiNewton(fg_im1c, zeros(n_c,1), lcpOpts);
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
                
                lcpOpts.b = bHi; lcpOpts.A = @(x) ALo*x;
                fgLo = @(x, Ax) quadraticLoss(x, ALo, bHi, Ax);
                xLo = proxQuasiNewton(fgLo, zeros(n,1),lcpOpts);
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
absErrFile = fullfile(basedir,'..','docs','fig', [prefix '_boundsHold.png']);
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
absErrFile = fullfile(basedir,'..','docs','fig', [prefix '_absErr.png']);
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
timePerHiFile = fullfile(basedir,'..','docs','fig', [prefix '_preCond.png']);
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
timePerHiFile = fullfile(basedir,'..','docs','fig', [prefix '_timePerHi.png']);
disp(['Saving to ' timePerHiFile])
saveas(f4, timePerHiFile);
end % plotDenseMats
function [hstar, zstar] = getC_A(A, debug)
if ~exist('debug','var') || isempty(debug)
    debug = false;
end
n = size(A,1);
[V,D] = eig(A);
[~,ixmin] = min(diag(D));
L = max(diag(D));
z0 = V(:,ixmin);
maxPGDIter = 1000;
lineSearchBudget = 100;
c1 =1e-4;
SGNS = [-1,1];
Hhat = zeros(n,2);
Zhat = cell(n,2);
for i = 1:n
    for ii = 1:2
        sigma = SGNS(ii);
        zkm1 = sigma*sign(z0(i))*z0 / abs(z0(i));
        zkm1 = max(-1,min(1, zkm1));
        h = zkm1.*(A*zkm1);
        [hkm1,j] = max(h);
        t = 1/L;
        for k = 1:maxPGDIter
            if debug; fprintf('- k=%d, j=%d, h=%.4g\n', k,j, hkm1); end;
            e = zeros(n,1);
            e(j) = 1;
            grad_h = (e*(e'*A))*zkm1 + ((A'*e)*e')*zkm1; % (sub)gradient of h
            % backtracing
            t = min(1/L,max(1,4*t));
            for l = 1:lineSearchBudget 
                zk = zkm1 - t * grad_h; % gradient descent step
                zk = sigma*sign(zk(i))*zk / abs(zk(i)); % project onto equality constraint
                zk = max(-1,min(1, zk)); % project onto inequality constraints
                h = max(zk.*(A*zk));
                p = zk - zkm1;
                if h < hkm1 + c1*t*dot(grad_h, p) % sufficient decrease
                    break
                end
                t = t / 2;
            end
            if l == lineSearchBudget
                if debug; fprintf('--- WARNING: Linesearch did not converge, l =%d\n', lineSearchBudget); end;
            end
            absErr = norm(zk - zkm1);
            relErr = absErr / norm(zkm1);
            if debug; fprintf('-- absErr=%.4g, relErr=%.4g\n', absErr, relErr); end;
            if absErr < 1e-8 
                if debug; fprintf('-- absErr < tol satisfied\n'); end;
                break
            elseif relErr < 1e-8
                if debug; fprintf('-- relErr < tol satisfied\n'); end;
                break
            end
            zkm1 = zk;
            h = zkm1.*(A*zkm1);
            [hkm1,j] = max(h);
        end
        if k == maxPGDIter
            if debug; fprintf('--- WARNING: PGD did not converge, k =%d\n', maxPGDIter); end;
        end
        Hhat(i,ii) = max(zk.*(A*zk));
        if debug; fprintf('hstar_%d^(%d) = %.4g\n', i, sigma,Hhat(i,ii)); end;
        Zhat{i,ii} = zk;
    end
end
Hhat = Hhat(:);
Zhat = Zhat(:);
[hstar, i] = min(Hhat);
if nargout > 1
    zstar = Zhat{i};
end
end % getC_A
%% getC_A_cvx
% AHi = (AHi + AHi') / 2;
% sqrtA = chol(AHi);
% cA = zeros(n,2);
% for m = 1:n
%     for mm = 1:2
%         if mm == 1
%             sigma = -1;
%         else
%             sigma = 1;
%         end
%         cvx_begin quiet
%             variable z(n)
%             y = sqrtA*z;
%             minimize max(y .* y)
%             subject to
%                 z <= 1
%                 -1 <= z
%                 z(m) == sigma
%         cvx_end
%         cA(m,mm) = cvx_optval;
%     end
% end
% cA = min(cA(:));



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
% timingFile = fullfile(basedir,'..','docs','fig', [prefix '_runTime.png']);
% disp(['Saving to ' timingFile])
% saveas(f1, timingFile);
% Relative Error