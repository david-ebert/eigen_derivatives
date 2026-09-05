% shared setup
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
n  = 2;
A0 = [1 0; 0 2];
A1 = [1 -1; -1 -1];
dA = {A1, zeros(n)};
M0 = eye(n);
dM = {zeros(n), zeros(n)};

[V, D] = eig(A0, "vector");
l0 = D(1);
u0 = V(:,1);
if u0(1) < 0
    u0 = -u0;   % fix the sign so the reference values are deterministic
end

%% eigenvalue derivatives match perturbation theory
[dL, ~] = eig_der(l0, u0, A0, dA);
assert(abs(dL{1} - 1) < 1e-12)
assert(abs(dL{2} + 2) < 1e-12)

%% second eigenvector derivative matches the normalised branch
[~, du] = eig_der(l0, u0, A0, dA);
assert(norm(du{2} - [-1; 4]) < 1e-10)

%% all three calling forms agree
[dL4, du4] = eig_der(l0, u0, A0, dA);
[dL5, du5] = eig_der(l0, u0, A0, dA, M0);
[dL6, du6] = eig_der(l0, u0, A0, dA, M0, dM);
for k = 1:numel(dL4)
    assert(norm(dL4{k} - dL5{k}) < 1e-12)
    assert(norm(dL4{k} - dL6{k}) < 1e-12)
    assert(norm(du4{k} - du5{k}) < 1e-12)
    assert(norm(du4{k} - du6{k}) < 1e-12)
end