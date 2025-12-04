f = figure();
metric = 'matVec';
hold on
edges = [1:1:14 15:5:50, opts.max_iter];
Ns = [];
MC = size(results,1);
numAlgo = size(results,2);
for ixAlgo = 1:numAlgo
    iters = zeros(MC,1);
    matVecs = zeros(MC,1);
    for mcIx = 1:numel(mcGood)
        mc = mcGood(mcIx);
        errHist = results(mc,ixAlgo).errHist;
        errHist(1,2) = NaN;
        try
            iters(mc) = find(errHist(:,1) < tol | ...
                errHist(:,2) < tol, 1, 'first');
            matVecs(mc) = errHist(iters(mc),3);
        catch 
             iters(mc) = opts.max_iter;
             matVecs(mc) = opts.max_iter;
        end
    end
    switch metric
        case 'iters'
            [N,edges] = histcounts(iters, edges);
            xaxis_name = 'Number of Iterations';
        case 'matVec'
            [N,edges] = histcounts(matVecs, edges);
            xaxis_name = 'Number of matVecs';
        otherwise
            error([metric ' is not a recognized metric'])
    end
    Ns = [Ns;N]; %#ok<AGROW>
end
bar(Ns');
xlabel(xaxis_name)
ylabel('Occurance');
title('Comparing LCP Solvers')
subtitle([...
    %'n \in [' num2str(minSz) ',' num2str(maxSz) '], 
    'reltol_{kkt} = ' num2str(tol) ', abstol_{kkt} = ' num2str(tol)] )
xticks(1:length(edges))
xticklabels([string(edges(1:14)), (string(edges(15:end-1)) + "-" +string(edges(16:end)))])
legend({results(mcGood(1),:).name})
fname = ['barChart_' prefix];
saveas(f, ['/Users/niru8088/scratch/Spheroidal3D-collisions/docs/fig/' fname '.pdf'])