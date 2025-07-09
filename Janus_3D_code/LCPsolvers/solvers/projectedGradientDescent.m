function [ x, info] = projectedGradientDescent(fg, x0, opts )
% May 2018 Wen Yan
% Edited 2025 Nic Rummel  
opts = defaultOpts(opts, numel(x0));

% Just a list of human readable text strings to convert the flag return
% code into something readable by writing msg(flag) onto the screen.
msg =  {'preprocessing';  % info.flag =  1
    'iterating';      % info.flag =  2
    'relative';       % info.flag =  3
    'absolute';       % info.flag =  4
    'stagnation';     % info.flag =  5
    'local minima';   % info.flag =  6
    'nondescent';     % info.flag =  7
    'maxlimit'        % info.flag =  8
    };

info.flag =  1;


max_iter = opts.max_iter;
tol_rel = opts.tol_rel;
tol_abs = opts.tol_abs; 
    
%--- Make sure all values are valid ---------------------------------------
max_iter = max(max_iter,1);
tol_rel  = max(tol_rel,0);
tol_abs  = max(tol_abs,0);
x0       = max(0,x0);


%--- Setup values need while iterating ------------------------------------
info = struct( ...
    'kkt', [], ...
    'iter',[],...
    'flag', [],...
    'msg',[] ...
); 
if ~isempty(opts.errFcn)
    if ishandle(opts.errFcn)
            info.errHist = zeros(max_iter+1,1);
            info.errHist(1) = opts.errFcn(x0);
    elseif iscell(opts.errFcn)
        info.errHist = zeros(max_iter+1,numel(opts.errFcn));
        for i = 1:numel(opts.errFcn)
            fcn = opts.errFcn{i};
            info.errHist(1,i) = fcn(x0);
        end
    end
end

err     = Inf;         % Current error measure
x       = x0;          % Current iterate
iter    = 1;           % Iteration counter

flag    = 2;

% first step, plain GD
[~,y] = fg(x);
[~,ytmp] = fg(y);
% NIC: is this a bug?
stepsize = (y'*y) / (y'*ytmp);

while (iter <= max_iter )
    xold = x;
    yold = y;
    x = xold - stepsize*yold; % new x
    x(x<0)=0; % projection
    
    %--- Test all stopping criteria used ------------------------------------
    phi     = min(y,x); 
    old_err = err;
    err     = 0.5*(phi'*phi);      % Natural merit function / kkt conditions
    
    if ~isempty(opts.errFcn)
        if ishandle(opts.errFcn)
            info.errHist(iter+1) = opts.errFcn(x);
        elseif iscell(opts.errFcn)
            for i = 1:numel(opts.errFcn)
                fcn = opts.errFcn{i};
                info.errHist(iter+1,i) = fcn(x);
            end
        end
    end

    if (abs(err - old_err) / abs(old_err)) < tol_rel  % Relative stopping criteria
        info.flag =  3;
        break;
    end
    if err < tol_abs   % Absolute stopping criteria
        info.flag =  4;
        break;
    end
    
    % Test if the search direction is smaller than numerical precision. That is if it is too close to zero.
    if stepsize < 10*eps
        info.flag =  5;
        % Rather than just giving up we may just use the gradient direction
        % instead. However, I am lazy here!
        break;
    end
    
    [~,y] = fg(x);
    
    % calculate B-B step size
    dx=x-xold;
    dy=y-yold;
    % BB1
    stepsize=(dx'*dx)/(dx'*dy);
    % BB2
    % stepsize=(dx'*dy)/(dy'*dy);
    % Increment the number of iterations
    iter=iter+1;
end

if iter>=max_iter
    info.flag =  8;
    info.iter = iter - 1;
end

info.kkt = err;
info.iter = iter;
info.msg =  msg{info.flag};
if ~isempty(opts.errFcn)
    info.errHist = info.errHist(1:iter+1,:);
end

end

