function testLCPWrappers()
try %#ok<TRYNC>
    rng('default')
end
rng(2);
mfilePath = mfilename('fullpath');
if contains(mfilePath,'LiveEditorEvaluationHelper')
    mfilePath = matlab.desktop.editor.getActiveFilename;
end
[dirname, ~,~] = fileparts(mfilePath);
addpath(genpath(fileparts(dirname)))
fname = 'amphi_5x5x5';
load([fname '.mat'], ...
    'A_list', 'b_list');
opts = struct( ...
    'max_iter',100, ...
    'tol_rel',1e-12, ...
    'tol_abs',1e-12, ... 
    'stepSizeRule', 'bb', ...
    'r', 20, ...
    'qnUpdate', 'bfgs' ...
);
MC = length(A_list);

algoNames = {'CVX', 'PGD', 'L-BFGS-B', 'Projected QuasiNewton (BFGS)', ...
    'Proximal QuasiNewton (BFGS)',...
    % 'Proximal zeroSR1', 'Proximal QuasiNewton (SR1)'...
    };
algoHndls = {@callCVX, @projectedGradientDescent, @L_BFGS_B, @projectedQuasiNewton, ...
    @proxQuasiNewton};
numAlgo = numel(algoNames);
assert(numel(algoNames) == numel(algoHndls));
results = repmat(...
    struct( ...
        'algo', '',...
        'time',   zeros(MC,1), ...
        'iters', zeros(MC,1), ...
        'kkt', zeros(MC,1), ...
        'matVecs', zeros(MC,1), ...
        'errHist', {cell(MC,1)} ...
    ), [1,numAlgo] ...
);
mcGood = [];
for mc = 50:MC
    disp(['mc = ' num2str(mc)])
    A = A_list{mc};
    Acnt = @(x) Acounter(x,A);
    b = b_list{mc};
    n = size(A,2);
    x0 = zeros(n,1);
    fg = @(x) objGrad(x, Acnt, b);
    opts.errFcn = {
        @(x) abs_kkt(x, A, b), 
        @(x) rel_kkt(x, A, b)
    };
    mcGoodFlag = true;
    for ixAlgo = 1:numAlgo
        name = algoNames{ixAlgo};
        algo = algoHndls{ixAlgo};
        tic
        % try 
            [x, info] = algo(fg, x0, opts);
        % catch
        %     if strcmpi(name, 'cvx')
        %         mcGoodFlag = false;
        %         break
        %     end
        % end
        results(ixAlgo).name = name;
        results(ixAlgo).time(mc) = toc();
        results(ixAlgo).iters(mc) = info.iter;
        results(ixAlgo).kkt(mc) = info.kkt;
        results(ixAlgo).matVecs(mc) = Acnt('reset');
        results(ixAlgo).errHist{mc} = info.errHist;
        if strcmpi(name, 'cvx')
            results(ixAlgo).matVecs(mc) = NaN;
            opts.errFcn{end+1} = @(xprime) rel_iter(xprime,x); 
            opts.errFcn{end+1} = @(xprime) abs_iter(xprime,x); 
        end
        for ixMetric = 1:numel(opts.errFcn)
            hndl = opts.errFcn{ixMetric};
            try %#ok<TRYNC>
                hndl('reset');
            end
        end
    end
    if mcGoodFlag
        mcGood = [mcGood mc];%#ok<AGROW>
        if results(2).iters(mc) < results(end).iters(mc)
            figure() 
            for ixAlgo = 2:numAlgo
                name = algoNames{ixAlgo};
                
                absKKT = results(ixAlgo).errHist{mc}(:,1);
                iter = numel(absKKT)-1;
                semilogy(0:iter, absKKT + 1e-13, 'LineWidth', 5)
                  hold on
            end
            legend(algoNames{2:end})
            xlabel('iterations')
            ylabel('kkt')
            title(['Iterations for MC ' num2str(mc)])
        end
    end
end
%%
fprintf('Algo                         | Time      | matVec | Iter | kkt\n')
for ixAlgo = 1:numAlgo
    name = algoNames{ixAlgo};
    while length(name) < length('Projected QuasiNewton (BFGS)')
        name = [name ' ']; %#ok<AGROW>
    end
    time = results(ixAlgo).time(mcGood);
    iters = results(ixAlgo).iters(mcGood);
    matVec = results(ixAlgo).matVecs(mcGood);
    kkt = results(ixAlgo).kkt(mcGood);
    fprintf('%s | %.1e s | %.2g\t| %.2g | %.2g\n', ...
        name, mean(time), mean(matVec), mean(iters), mean(kkt));
end
%% 
% f = figure();
% hold on
% edges = [1:1:9 10:5:50, 100];
% [N,edges] = histcounts(iter_bbpgd(mcGood),edges);
% [N1,edges] = histcounts(iter_lbfgsb(mcGood),edges);
% [N2,edges] = histcounts(iter_plbfgs(mcGood),edges);
% [N3,edges] = histcounts(iter_zerosr1(mcGood),edges);
% [N4,edges] = histcounts(iter_pq(mcGood),edges);
% bar([N;N1;N2;N3;N4]');
% xlabel('Number of Iterations');
% ylabel('Occurance');
% tt = join(split(fname, '_'), ' '); 
% tt = tt{1};
% title(['Comparing LCP solvers for ' tt])
% xticks(1:length(edges))
% xticklabels([string(edges(1:9)), (string(edges(10:end-1)) + "-" +string(edges(11:end)))])
% legend({'BB-PGD', 'L-BFGS-B', 'P-L-BFGS' 'zerosr1', 'proxQuasNewton'})
% saveas(f, ['/Users/niru8088/scratch/Spheroidal3D-collisions/docs/' fname '.pdf'])
end

function Ax= Acounter(x,A)
persistent matVecCnt 

if isempty(matVecCnt)
    matVecCnt = 0;
end

if ischar(x) 
    if strcmpi(x, 'matVecCnt')
        Ax = matVecCnt; 
        return 
    elseif strcmpi(x, 'reset')
        Ax = matVecCnt; 
        matVecCnt = 0;
        return 
    else
        assert(false, ['Option: ' x ' not recognized'])
    end
end

Ax = A*x;
matVecCnt = matVecCnt + 1;
end

function [f,g] = objGrad(x,A,b)

Ax = A(x);
f = 1/2 * dot(x, Ax) + dot(x, b);
if nargout == 1
    g = [];
    return 
end
g = Ax + b;
end


function [x, info] = callCVX(fg, x0, opts)

n = length(x0);
cvx_begin quiet
        variable x(n)
        minimize fg(x) 
        subject to 
        0 <= x
cvx_end 

info.iter = NaN;
info.kkt = NaN;
info.errHist = NaN;

end

% function [f,g] = objGradWithHist(fg, x, errFcn)
%     persistent errHist
% 
%     if ischar(x)
%         if strcmpi(x, 'errHist')
%             f = errHist;
%             g = [];
%             return 
%         elseif strcmpi(x, 'reset')
%             f = errHist;
%             g = [];
%             if ishandle(errFcn)
%                 errHist = zeros(0,1);
%             elseif iscell(errFcn)
%                 errHist = zeros(0,numel(errFcn));
%             end
%             return 
%         end
%     end
%     [f,g] = fg(x);
%     if ~isempty(errFcn)
%         if isempty(errHist)
%             if ishandle(errFcn)
%                 errHist = zeros(0,1);
%             elseif iscell(errFcn)
%                 errHist = zeros(0,numel(errFcn));
%             end
%         end
% 
%         if ishandle(errFcn)
%             errHist(end+1) = errFcn(x);
%         elseif iscell(errFcn)
%             for i = 1:numel(errFcn)
%                 fcn = errFcn{i};
%                 xtmp = x;
%                 errHist(end+1,i) = fcn(xtmp); %#ok<AGROW>
%             end
%         end
%     end
% end

function e = abs_kkt(x, A, b)
    g = A*x + b;
    phi = min(x,g);
    e = 1/2*dot(phi, phi);
end % abs_kkt

function diff = rel_kkt(x, A, b)
    persistent olde 
    if ischar(x) && strcmpi(x,'reset')
        olde = [];
        diff = NaN;
        return
    end
    g = A*x + b;
    phi = min(x,g);
    e = 1/2*dot(phi, phi);
    if isempty(olde) 
        diff = NaN;
    else 
        diff = abs(e - olde) / abs(e);
    end
    olde = e;
end % rel_kkt

function e = abs_iter(x, xstar)
e = norm(x - xstar);
end % abs_iter

function e = rel_iter(x, xstar)
e = norm(x-xstar) / norm(xstar);
end % rel_iter