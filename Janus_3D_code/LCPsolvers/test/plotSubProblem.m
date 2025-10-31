function plotSubProblem(results, A, vals, plotFlag)
if ~exist('vals','var') || isempty(vals)
    vals = eig(A);
end 
if ~exist('plotFlag','var') || isempty(plotFlag)
    plotFlag = true;
end 
if ~plotFlag
    return 
end
linespec = {"-o", "-s","-*","-diamond","-^","-+" };
plotlyjs_colors = {"#1f77b4";  % muted blue
    "#ff7f0e";  % safety orange
    "#2ca02c";  % cooked asparagus green
    "#d62728";  % brick red
    "#9467bd";  % muted purple
    "#8c564b";  % chestnut brown
    "#e377c2";  % raspberry yogurt pink
    "#7f7f7f";  % middle gray
    "#bcbd22";  % curry yellow-green
    "#17becf"; % blue-teal
    };  

objVal = cell(numel(results),1);
for ixAlgo = 1:numel(results)
    iterHist = results(ixAlgo).iterHist;
    errHist = results(ixAlgo).errHist;
    numIter = size(errHist,1)-1;
    try
        objVal{ixAlgo} = arrayfun(@(i) fg(iterHist(i,:)',[]), 1:numIter+1);
    catch
        objVal{ixAlgo} = errHist(:,end);
    end
end
% minObjVal = min(cellfun(@min, objVal));
figure()
for ixAlgo = 1:numel(results)
    % name = algoNames{ixAlgo};
    errHist = results(ixAlgo).errHist;
    % numIter = size(errHist,1)-1;
    kk = mod(ixAlgo-1, numel(linespec)) + 1;
    jj = mod(ixAlgo-1, numel(plotlyjs_colors)) + 1;
    %% objective value top left
    % subplot(3,1,1)
    % hold on
    % semilogy(0:numIter, objVal{ixAlgo} - minObjVal + 1e-12, linespec{kk}, 'LineWidth',4,'MarkerSize',5, 'Color', plotlyjs_colors{kk}) % abs(objVal - cvxObjVal) / abs(cvxObjVal))
    % xlabel('Iteration','FontSize', 20)
    % ylabel('$f(x^{(k)})$','interpreter', 'latex', 'FontSize', 25)
    % set(gca,'YScale','log')
    %% Iterate Error Top Middle
    % subplot(3,2,3)
    % hold on 
    % semilogy(0:numIter, errHist(:,4), linespec{kk},'LineWidth',4,'MarkerSize',10,'Color', plotlyjs_colors{kk});
    % ylabel('$\frac{\|x^{(k)} - x^*\|}{\|x^*\|}$', 'interpreter', 'latex', 'FontSize', 30)
    % xlabel('Iteration','FontSize', 20)
    % set(gca,'YScale','log')
    %% KKT condition bottom
    subplot(1,3,[1 2])
    hold on
    % semilogy(0:numIter, errHist(:,1), linespec{kk}, 'LineWidth',4, 'MarkerSize',10, 'Color', plotlyjs_colors{kk});

    semilogy(errHist(:,3), errHist(:,1), linespec{kk}, 'LineWidth',4, 'MarkerSize',10, 'Color', plotlyjs_colors{jj});
    xlabel('MVPs','FontSize', 20)
    ylabel('$\min(Ax^{(k)}+b, x^{(k)})$','interpreter', 'latex','FontSize', 25)
    set(gca,'YScale','log')
end
%% Legend Top Left 
% subplot(1,2,1)
% legend({results.name},'FontSize', 20,'Location','northeastoutside')
legend({results.name},'FontSize', 20,'Location','northeast')
%% Eigen values Middle Right
subplot(1,3,3)
semilogy(sort(vals,1,'descend'), '-o', 'Color', '#808080', ...
    'LineWidth',4, 'MarkerSize',10)
ylabel('Eigen Values of $A$', 'interpreter', 'latex','FontSize', 25)
%% Large Title
sgtitle({'Exmample Problem', ...
    sprintf('n = %d, \\kappa(A) = %.2g', ...
    size(A,1),  cond(A))},'FontSize', 30)
end