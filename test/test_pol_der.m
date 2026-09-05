% shared setup
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
angle_x = pi/4;
dx = cos(angle_x);
dy = sin(angle_x);
dxK = [1 0; 0 -1];
dyK = [0 1; 1 0];

%% cross_friswell: split in first order, only the last order stays undetermined
n = 2;
K_0 = 2*eye(n);
dK = {dx*dxK + dy*dyK, zeros(n), zeros(n), zeros(n)};
[V, D] = eig(K_0, "vector");
[dL, du] = eig_der(D(1), V, K_0, dK);
[P0, k_pol] = pol(dL);
[dP, dl] = pol_der(dL, P0, k_pol, V, du);

assert(norm(dl{1} - [-1; 1]) < 1e-10)
assert(norm(dl{2}) < 1e-10)
assert(norm(dl{3}) < 1e-10)
assert(numel(dP) == numel(dK))
assert(all(isfinite(dP{1}(:))))
assert(all(isfinite(dP{2}(:))))
assert(all(isfinite(dP{3}(:))))
assert(all(isnan(dP{4}(:))))

%% deflect: split in second order, the last two orders stay undetermined
n = 2;
K_0 = 2*eye(n);
dK = {zeros(n), dx^2*dxK + dy^2*dyK, zeros(n), zeros(n)};
[V, D] = eig(K_0, "vector");
[dL, du] = eig_der(D(1), V, K_0, dK);
[P0, k_pol] = pol(dL);
[dP, dl] = pol_der(dL, P0, k_pol, V, du);

assert(isequal(k_pol, 2*ones(n)))
assert(norm(dl{1}) < 1e-10)
assert(norm(dl{2} - [-1/sqrt(2); 1/sqrt(2)]) < 1e-10)
assert(all(isfinite(dP{1}(:))))
assert(all(isfinite(dP{2}(:))))
assert(all(isnan(dP{3}(:))))
assert(all(isnan(dP{4}(:))))

%% du_pol with a non-constant polarization follows dP
n = 2;
K_0 = 2*eye(n);
dK = {dx*dxK + dy*dyK, zeros(n), zeros(n), zeros(n)};
[V, D] = eig(K_0, "vector");
[dL, du] = eig_der(D(1), V, K_0, dK);
[P0, k_pol] = pol(dL);
[dP, ~] = pol_der(dL, P0, k_pol, V, du);
[du_p, u0_p] = du_pol(V, du, P0, dP);

assert(numel(du_p) == numel(dP))
assert(norm(u0_p - V*P0) < 1e-12)
assert(all(isfinite(du_p{1}(:))))
assert(all(isnan(du_p{4}(:))))

%% du_pol with a constant polarization is a plain right multiplication
n = 2;
K_0 = 2*eye(n);
dK = {dx*dxK + dy*dyK, zeros(n), zeros(n), zeros(n)};
[V, D] = eig(K_0, "vector");
[~, du] = eig_der(D(1), V, K_0, dK);
P0 = eye(n);
[du_p, u0_p] = du_pol(V, du, P0);

assert(numel(du_p) == numel(du))
assert(norm(u0_p - V*P0) < 1e-12)
for i = 1:numel(du)
    assert(norm(du_p{i} - du{i}*P0) < 1e-12)
end

%% pol returns a full matrix of Inf when the derivatives run out
n = 2;
previous = warning('off', 'eig_deriv:pol_unseparable');
[P0, k_pol] = pol({eye(n)});
warning(previous);

assert(isequal(size(k_pol), [n n]))
assert(all(isinf(k_pol(:))))
assert(isequal(size(P0), [n n]))
