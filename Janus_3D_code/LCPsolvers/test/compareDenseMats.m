%% set paths
[dirname, ~] = setPaths();
%% set srcFilie
srcDir = fullfile(dirname, '../data/amphiLCPs.n_2.p_8.cDist_2.3');
srcFile = fullfile(srcDir,'allMats.mat');
%%
load(srcFile, 'ps', 'gmresTols', 'out');
%%
K = numel(ps);
L = numel(gmresTols);
N = numel(out);
N = 69;
relErr = zeros(K, L, N);
eigVals = cell(K, L, N);
for n = 1:N
    A_hifi = out{n}{1,1};
    A = 1/2*(A_hifi + A_hifi');
    for k = 1:K
        for l = 1:L
            Ankl = out{n}{k,l};
            try
                relErr(k,l,n) = norm(A - Ankl) / norm(A);
            catch
                relErr(k,l,n) = nan;
            end
            try
                eigVals{k,l,n} = sort(eig(Ankl), 'descend');
            catch 
                eigVals{k,l,n} = [];
            end
        end
    end
end
%%
mean_rel = zeros(K,L);
std_rel = zeros(K,L);
mean_eig = cell(K,L);
std_eig = cell(K,L);
for k = 1:K
    for l = 1:L
        % relErr
        mean_rel(k,l) = nanmean(relErr(k,l,:));
        std_rel(k,l) = nanstd(relErr(k,l,:));
        % eigVals
        M = max(cellfun(@numel, eigVals(k,l,:)));
        mean_eig{k,l} = zeros(M,1);
        std_eig{k,l} = zeros(M,1);
        cnts = zeros(M,1);
        for n = 1:N
            m = numel(eigVals{k,l,n});
            mean_eig{k,l}(1:m) = mean_eig{k,l}(1:m) + real(eigVals{k,l,n});
            cnts(1:m) = cnts(1:m) +1;
        end
        mean_eig{k,l} = mean_eig{k,l} ./ cnts;
        for n = 1:N
            m = numel(eigVals{k,l, n});
            std_eig{k,l}(1:m) = std_eig{k,l}(1:m) + (mean_eig{k,l}(1:m) - eigVals{k,l, n}).^2;
        end
        std_eig{k,l} = real(sqrt(std_eig{k,l} ./ (cnts - 1)));
    end
end
%% Relative error 
f1 = figure;
for k = 1:K
    for l = 1:L
        subplot(L, K, K*(l-1) + k);
        edges = 10.^(-6:.5:2);
        [counts,edges] = histcounts(relErr(k,l,:),edges);
        g = histogram('BinEdges',edges,'BinCounts',counts);
        set(gca, "Xscale", "log")
        xticks(10.^(-6:1:2))
        title({ ...
            sprintf('p = %d, tol = %.0e', ps(k), gmresTols(l)), ...
            sprintf('mean = %.2e',mean_rel(k,l)),...
            sprintf('std = %.2e', mean_rel(k,l))}...
        );
    end
end
[path,~,~] = fileparts(srcFile);
[~,name,~] = fileparts(path);
prts = split(name, '_');
name = strjoin(prts, ' = ');
prts = split(name, '.');
name = strjoin(prts, ', ');
sgtitle({'Relative Eror',name})
set(f1, 'Position',  [0, 0, 1000, 1200])
saveas(f1, fullfile(srcDir, [name '_relErr.pdf']));
%%
f2 = figure;
for k = 1:K
    for l = 1:L
        subplot(L, K, K*(l-1) + k);
        hold on
        M = numel(mean_eig{k,l});
        errorbar(1:M, mean_eig{k,l}, std_eig{k,l});
        title(sprintf('p = %d, tol = %.0e', ps(k), gmresTols(l)));
    end
end
[path,~,~] = fileparts(srcFile);
[~,name,~] = fileparts(path);
prts = split(name, '_');
name = strjoin(prts, ' = ');
prts = split(name, '.');
name = strjoin(prts, ', ');
sgtitle({'Eigen Values',name})
set(f2, 'Position',  [0, 0, 1000, 1200])
saveas(f2, fullfile(srcDir, [name '_eigVals.pdf']));
