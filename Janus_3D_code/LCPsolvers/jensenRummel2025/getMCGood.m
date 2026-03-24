function mcGood = getMCGood(resFile, ps, tols)
persistent cache 
if isempty(cache)
    cache = containers.Map('KeyType','char', 'ValueType', 'any');
end
disp('Selecting subset of time steps where no errors occurred')
if cache.isKey(resFile)
    disp('- Loading from precomputed cache')
    mcGood = cache(resFile);
    return 
elseif isfile([resFile  '.mcGood.mat'])
    load([resFile  '.mcGood.mat'], 'mcGood')
    cache(resFile) = mcGood;
    return 
end
disp('- Computing from scratch')
res = load(resFile);
Nt = size(res.A,1);
mcGood = true(Nt, 1);
if ~exist('ps', 'var') || isempty(ps)
    ps = res.ps;
    tols = repmat(min(res.tols), [numel(ps), 1]);
end
for i = 1:Nt
    disp(['-- i = ' num2str(i) '/' num2str(Nt)])
    for l = 1:numel(ps)
        p = ps(l);
        tol = tols(l);
        j = find(res.ps == p,1,'first');
        k = find(abs(res.tols - tol)/abs(tol) < 1e-8,1,"first");
        A_ = res.A{i,j,k};
        if isempty(A_) || norm(A_) > 10
            mcGood(i) = false;
            break
        end
        try 
            n = size(A_,1);
            callCVX(zeros(n,1), A_, res.b{i});
        catch 
            mcGood(i) = false;
            break
        end
    end
end
mcGood = find(mcGood);
save(resFile + '.mcGood.mat', 'mcGood')
cache(resFile) = mcGood;