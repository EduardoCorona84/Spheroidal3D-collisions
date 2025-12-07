function [h0, rho, S, Y, opts] = updateQNMemory(s, y, opts)
persistent looseMem
if isempty(looseMem)
    looseMem = false;
end
h0 = 1; 
n = opts.n;
if isinf(opts.qn.m)
    % Corresponds to full memory
    m = n;
else
    m = opts.qn.m;
end
rho = opts.qn.rho;
S = opts.qn.S;
Y = opts.qn.Y;
if isempty(s) 
    r = find(arrayfun(@(i) any(S(:,i) ~= 0), 1:m),1,'last');
    rho = rho(1:r);
    S = S(:,1:r);
    Y = Y(:,1:r);
    if ~isempty(r)
        h0 = geth0(S(:,r), Y(:, r));
    end 
    looseMem = false;
    return 
end 

%%  Loose memory if the previous update was skipped
% or if the memory is full
if (looseMem && any(S(:)~=0)) || all(arrayfun(@(i) any(S(:,i) ~= 0), 1:m))
    r = find(arrayfun(@(i) any(S(:,i) ~= 0), 1:m),1,'last');
    rho(1:r-1) = rho(2:r);
    S(:,1:r-1) = S(:,2:r);
    Y(:,1:r-1) = Y(:,2:r);
    rho(r) = 0;
    S(:,r) = 0;
    Y(:,r) = 0;
end
%% Set the index to 
r = find(arrayfun(@(i) all(S(:,i) == 0), 1:m),1,'first');
assert(~isempty(r), 'r should not be empty at this point')
% Curvature check on secant conditions
if dot(s,y) >= 1e-8
    S(:, r) = s;
    Y(:, r) = y;
    rho(r) = 1/dot(s,y);
    looseMem = false;
else
    % during failure use the previous information, 
    % and set flag to loose memory
    looseMem = true;
    r = r - 1;
end
if r <= 0
    rho = [];
    S = [];
    Y = [];
    return 
end

% Set h0
h0 = geth0(S(:,r), Y(:, r));
opts.qn.rho = rho;
opts.qn.S = S;
opts.qn.Y = Y;

rho = rho(1:r);
S = S(:, 1:r);
Y = Y(:, 1:r);
end % updateQNMemory

function h0 = geth0(s,y)
tau_bb2 = dot(s,y) / dot(y,y); % dot(s,y) / norm(y,2)^2
%% From Stephens ProxQN paper
gamma = 0.8;
tau_min = 1e-14;
tau_max = Inf;
tau_bb2 = min(max(tau_bb2, tau_min), tau_max);
if tau_bb2 == tau_min
    warning('Convexity of cost function is stagnating'); 
end
h0 = gamma * tau_bb2;
end