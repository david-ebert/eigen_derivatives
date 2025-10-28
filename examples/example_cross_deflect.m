clear; restoredefaultpath
%%

n_inner = 3;
n_ev    = 3;
degenerate = true;

x0 = 0; y0 = 0;

K_0 = eye(n_inner)*2;
dxK = ...
    [ 1  ,  1  , 1 ;
      1  ,  1  , 1 ;
      1  ,  1  , 1  ];
dyK = ...
    [ 0  ,  0  , 0 ;
      0  ,  0  , 0 ;
      0  ,  0  , 0  ];
dx2K = ...
    [ 0  ,  0  , 0 ;
      0  ,  1  , 0 ;
      0  ,  0  , -1  ]*2;
dy2K = ...
    [ 0  ,  0  , 0 ;
      0  ,  0  , 1 ;
      0  ,  1  , 0  ]*2;
dx3K = ...
    [ 0  ,  0  , 0 ;
      0  ,  1  , 0 ;
      0  ,  0  , -1  ]*6;
dy3K = ...
    [ 0  ,  0  , 0 ;
      0  ,  0  , 1 ;
      0  ,  1  , 0  ]*6;

dx4K = zeros(n_inner);
dy4K = zeros(n_inner);

%%
demo