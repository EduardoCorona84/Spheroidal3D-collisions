function [h0, rho, S, Y, opts] = updateQNMemory(s, y, opts)
persistent looseMem
if isempty(looseMem)
    looseMem = false;
end
h0 = 1; 
rho = [];
S = [];
Y = [];
n = opts.n;
m = opts.qn.m;
rho_ = opts.qn.rho;
S_ = opts.qn.S;
Y_ = opts.qn.Y;
if isempty(s) 
    r = find(arrayfun(@(i) any(S_(:,i) ~= 0), 1:m),1,'last');
    if ~isempty(r)
        rho = rho_(1:r);
        S = S_(:,1:r);
        Y = Y_(:,1:r);
        h0 = geth0(S(:,r), Y(:, r));
    end 
    looseMem = false;
    return 
end 

%%  Loose memory if the previous update was skipped
% or if the memory is full
if (looseMem && any(S_(:)~=0)) || all(arrayfun(@(i) any(S_(:,i) ~= 0), 1:m))
    r = find(arrayfun(@(i) any(S_(:,i) ~= 0), 1:m),1,'last');
    rho_(1:r-1) = rho_(2:r);
    S_(:,1:r-1) = S_(:,2:r);
    Y_(:,1:r-1) = Y_(:,2:r);
    rho_(r) = 0;
    S_(:,r) = 0;
    Y_(:,r) = 0;
end
%% Set the index to 
r = find(arrayfun(@(i) all(S_(:,i) == 0), 1:m),1,'first');
assert(~isempty(r), 'r should not be empty at this point')
% Curvature check on secant conditions
if dot(s,y) >= 1e-8
    S_(:, r) = s;
    Y_(:, r) = y;
    rho_(r) = 1/dot(s,y);
    looseMem = false;
else
    % during failure use the previous information, 
    % and set flag to loose memory
    looseMem = true;
    r = r - 1;
end
if r <= 0
    return 
end
rho = rho_(1:r);
S = S_(:, 1:r);
Y = Y_(:, 1:r);
% Set h0
h0 = geth0(S(:,r), Y(:, r));

end % updateQNMemory

function h0 = geth0(s,y)
tau_bb2 = dot(s,y) / dot(y,y); % dot(s,y) / norm(y,2)^2
%% From Stephens ProxQN paper
gamma = 0.8;
tau_min = 1e-14;
tau_max = Inf;
tau_bb2 = clip(tau_bb2, tau_min, tau_max);
if tau_bb2 == tau_min
    warning('Convexity of cost function is stagnating'); 
end
h0 = gamma * tau_bb2;
end