f = figure();
abs_kkt = 1e-6; 
rel_kkt = 1e-6; 
metric = 'matVec';
hold on
edges = [1:1:14 15:5:50, 100];
Ns = [];
MC = length(results(1).iters);
for i = 2:length(results)
    iters = zeros(MC,1);
    matVecs = zeros(MC,1);
    for mc = 50:100
        errHist = results(i).errHist{mc};
        errHist(1,2) = NaN;
        iters(mc) = find(errHist(:,1) < abs_kkt | ...
            errHist(:,2) < rel_kkt, 1, 'first');
        matVecs(mc) = errHist(iters(mc),3);
    end
    switch metric
        case 'iters'
            [N,edges] = histcounts(iters, edges);
            xlabel('Number of Iterations');
        case 'matVec'
            [N,edges] = histcounts(matVecs, edges);
            xlabel('Number of matVecs');
        otherwise
            error([metric ' is not a recognized metric'])
    end
    Ns = [Ns;N]; %#ok<AGROW>
end
bar(Ns');


ylabel('Occurance');
tt = join(split(fname, '_'), ' '); 
tt = tt{1};
title(['Comparing LCP solvers for ' tt])
xticks(1:length(edges))
xticklabels([string(edges(1:14)), (string(edges(15:end-1)) + "-" +string(edges(16:end)))])
legend(algoNames{2:end})
% saveas(f, ['/Users/niru8088/scratch/Spheroidal3D-collisions/docs/' fname '.pdf'])