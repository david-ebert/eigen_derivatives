% degenerate eigenvalues (multiplicity 2),
% such that the trajectories deflect,
% first derivatives of eigenvalues coincide, second derivatives differ

clear; restoredefaultpath
%%

n_inner = 2;
n_ev    = 2;
degenerate = true;

x0 = 0; y0 = 0;

K_0 = eye(n_inner)*2;

dxK = zeros(n_inner);
dyK = zeros(n_inner);

dx2K = ...
    [ 1 , 0 ;
      0 ,-1 ];
dy2K = ...
    [ 0 , 1 ;
      1 , 0 ];

dx3K = zeros(n_inner);
dx4K = zeros(n_inner);
dy3K = zeros(n_inner);
dy4K = zeros(n_inner);

%%
demo