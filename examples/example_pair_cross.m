% degenerate eigenvalues (multiplicity 4),
% such that the eigenvalues form two deflecting pairs,
% which as pairs cross each other

clear; restoredefaultpath
%%

n_inner = 4;
n_ev    = 4;
degenerate = true;

x0 = 0; y0 = 0;

K_0 = eye(n_inner)*2;
dxK = ...
    [ 1 , 1 , 0 , 0 ;
      1 , 1 , 0 , 0 ;
      0 , 0 , 2 , 0 ;
      0 , 0 , 0 , 0 ];
dyK = ...
    [ 0 , 0 , 0 , 0;
      0 , 0 , 0 , 0;
      0 , 0 , 0 , 0;
      0 , 0 , 0 , 0; ];
dx2K = ...
    [ 1 , 1 , 0 , 0 ;
      1 , 1 , 0 , 0 ;
      0 , 0 , 1 , 1 ;
      0 , 0 , 1 , 1  ]*2;
dy2K = ...
    [ 0 , 1 , 0 , 0 ;
      1 , 0 , 0 , 0 ;
      0 , 0 , 0 , 1 ;
      0 , 0 , 1 , 0  ]*2;

dx3K = zeros(n_inner);
dx4K = zeros(n_inner);
dy3K = zeros(n_inner);
dy4K = zeros(n_inner);

%%
demo