function flow_struct = SpheroidalMS_EvalBackgroundFlow(parbd, background_flow)
% The background velocity is
%   u_inf(x) = U0 + A * x
% or equivalently written as
%   X * A.' + U0
%
% The associated stress is (setting the associated pressure to be 0)
%   sigma_{\infty} = mu * (A + A.')
% and the traction on the surface is sigma_{\infty} * n(x)

num_nodes = size(parbd.Xp, 1);
zero_vec = zeros(3*num_nodes, 1);

flow_struct = struct( ...
    'enabled', false, ...
    'velocity', zero_vec, ...
    'traction', zero_vec, ...
    'stress', zeros(3,3) ...
);

if nargin < 2 || isempty(background_flow) || ~background_flow.enabled
    return;
end

if isfield(parbd, 'mu')
    fluid_mu = parbd.mu;
else
    fluid_mu = 1;
end

U0 = reshape(background_flow.U0, 1, 3);
A = background_flow.A;
stress = fluid_mu * (A + A.');

flow_struct.enabled = true;
flow_struct.velocity = reshape((parbd.Xp * A.' + U0).', [], 1);
flow_struct.traction = reshape((parbd.Nrp * stress.').', [], 1);
flow_struct.stress = stress;
end
