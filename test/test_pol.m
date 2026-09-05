% shared setup
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
angle_x = pi/4;
dx = cos(angle_x);
dy = sin(angle_x);

%% cross_friswell: cones, eigenvalues split at first order
n = 2;
K_0 = 2*eye(n);
dxK = [1 0; 0 -1];
dyK = [0 1; 1 0];
dK = {dx*dxK + dy*dyK, zeros(n), zeros(n), zeros(n)};
[V, D] = eig(K_0, "vector");
[dL, ~] = eig_der(D(1), V, K_0, dK);
[~, k_pol] = pol(dL);
assert(isequal(k_pol, ones(n)))

%% deflect: first derivatives coincide, split at second order
n = 2;
K_0 = 2*eye(n);
dx2K = [1 0; 0 -1];
dy2K = [0 1; 1 0];
dK = {zeros(n), dx^2*dx2K + dy^2*dy2K, zeros(n), zeros(n)};
[V, D] = eig(K_0, "vector");
[dL, ~] = eig_der(D(1), V, K_0, dK);
[~, k_pol] = pol(dL);
assert(isequal(k_pol, 2*ones(n)))

%% cross_deflect: one eigenvalue splits at first order, the remaining pair at second
n = 3;
K_0 = 2*eye(n);
dxK  = ones(n);
dx2K = diag([0 1 -1])*2;
dy2K = [0 0 0; 0 0 1; 0 1 0]*2;
dx3K = diag([0 1 -1])*6;
dy3K = [0 0 0; 0 0 1; 0 1 0]*6;
dK = {dx*dxK, dx^2*dx2K + dy^2*dy2K, dx^3*dx3K + dy^3*dy3K, zeros(n)};
[V, D] = eig(K_0, "vector");
[dL, ~] = eig_der(D(1), V, K_0, dK);
[~, k_pol] = pol(dL);
assert(isequal(k_pol, [2 2 1; 2 2 1; 1 1 1]))

%% cross_cross: one eigenvalue splits at first order, the remaining pair at third
n = 3;
K_0 = 2*eye(n);
dxK  = ones(n);
dx3K = diag([0 1 -1])*2;
dy3K = [0 0 0; 0 0 1; 0 1 0]*2;
dx4K = diag([0 1 -1])*6;
dy4K = [0 0 0; 0 0 1; 0 1 0]*6;
dK = {dx*dxK, zeros(n), dx^3*dx3K + dy^3*dy3K, dx^4*dx4K + dy^4*dy4K};
[V, D] = eig(K_0, "vector");
[dL, ~] = eig_der(D(1), V, K_0, dK);
[~, k_pol] = pol(dL);
assert(isequal(k_pol, [3 3 1; 3 3 1; 1 1 1]))

%% pair_cross: two pairs split at second order, the pairs from each other at first
n = 4;
K_0 = 2*eye(n);
dxK  = [1 1 0 0; 1 1 0 0; 0 0 2 0; 0 0 0 0];
dx2K = [1 1 0 0; 1 1 0 0; 0 0 1 1; 0 0 1 1]*2;
dy2K = [0 1 0 0; 1 0 0 0; 0 0 0 1; 0 0 1 0]*2;
dK = {dx*dxK, dx^2*dx2K + dy^2*dy2K, zeros(n), zeros(n)};
[V, D] = eig(K_0, "vector");
[dL, ~] = eig_der(D(1), V, K_0, dK);
[~, k_pol] = pol(dL);
assert(isequal(k_pol, [2 2 1 1; 2 2 1 1; 1 1 2 2; 1 1 2 2]))

%% exhausted derivatives: only the degenerate blocks are marked
% two pairs that stay degenerate, but the pairs separate from each other
% at first order, so those entries must survive
dL = {diag([1 1 2 2])};
warning('off', 'eig_deriv:pol_unseparable');
[~, k_pol] = pol(dL);
warning('on', 'eig_deriv:pol_unseparable');
assert(isequal(k_pol, [Inf Inf 1 1; Inf Inf 1 1; 1 1 Inf Inf; 1 1 Inf Inf]))
