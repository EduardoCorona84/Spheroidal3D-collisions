function mcGood = getMCGood(resFile)
persistent cache 
if isempty(cache)
    cache = containers.Map('KeyType','char', 'ValueType', 'any');
end
disp('Selecting subset of time steps where no errors occurred')
if cache.isKey(resFile)
    disp('- Loading from precomputed cache')
    mcGood = cache(resFile);
    return 
end
disp('- Computing from scratch')
res = load(resFile);
Nt = size(res.A,1);
mcGood = true(Nt, 1);
ps = [8,6,4,3];
tols = [1e-8,1e-6,1e-6,1e-5];
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
cache(resFile) = mcGood;