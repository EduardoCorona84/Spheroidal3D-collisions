function plotSubProblem(results, A, plotFlag, ttlStr, vals)
if ~exist('plotFlag','var') || isempty(plotFlag)
    plotFlag = true;
end 
if ~exist('ttlStr','var') || isempty(ttlStr)
    ttlStr = 'Example Problem';
end 
if ~exist('vals','var') || isempty(vals)
    vals = eig(A);
end 
if ~plotFlag
    return 
end
markers = {"o", "s","*","diamond","^","+" };
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
gcf
clf
ps = 3:7;
tols = [1e-5, 1e-6, 1e-7, 1e-8];
dashOpts= {"-", "--", ":", "-."};
for ixAlgo = 1:numel(results)
    name = results(ixAlgo).name;
    errHist = results(ixAlgo).errHist;
    % numIter = size(errHist,1)-1;
    if contains(name, 'B_0')
        prts = split(name, 'p=');
        prts =split(prts{end},',');
        p = str2double(prts{1});
        jj = 2+ find(p == ps, 1,'first');
        prts =split(prts{end},'=');
        prts = split(prts{end},'$)');
        tol = str2double(prts{1});
        % dash = dashOpts{tol == tols}; 
        if contains(name, 'B(')
            kk = 3;
            dash = '-';
        elseif contains(name, '\%') 
            kk =4;
            dash = ':';
        elseif contains(name, 'x_0') 
            kk = 6;
            dash = '-.';
        else 
            kk =5;
            dash = '--';
        end
        
    else
        dash = '-';
        jj = mod(ixAlgo-1, numel(plotlyjs_colors)) + 1;
        kk = mod(ixAlgo-1, numel(markers)) + 1;
    end
    
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

    plot(errHist(:,3), errHist(:,1), 'LineStyle', dash, ...
        'Marker',markers{kk}, 'LineWidth',4, 'MarkerSize',10, ...
        'Color', plotlyjs_colors{jj});
    xlabel('MVPs','FontSize', 20)
    ylabel('$\min(Ax^{(k)}+b, x^{(k)})$','interpreter', 'latex','FontSize', 25)
    set(gca,'YScale','log')
end
%% Legend Top Left 
% subplot(1,2,1)
% legend({results.name},'FontSize', 20,'Location','northeastoutside')
% xlim([0,25])
legend({results.name},'FontSize', 20,'Location','northeast', 'interpreter', 'latex')
%% Eigen values Middle Right
subplot(1,3,3)
semilogy(sort(vals,1,'descend'), '-o', 'Color', '#808080', ...
    'LineWidth',4, 'MarkerSize',10)

ylabel('Eigen Values of $A$', 'interpreter', 'latex','FontSize', 25)
%% Large Title
sgtitle(ttlStr,'FontSize', 30, 'interpreter', 'latex')
title(sprintf('n = %d, \\kappa(A) = %.2g', size(A,1),  cond(A)),'FontSize', 20)
end