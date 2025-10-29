% non-degenerate eigenvalue
%
% See also EXAMPLE_CROSS_FRISWELL
%
% Friswell, M. I. (July 1, 1996). 
% "The Derivatives of Repeated Eigenvalues and Their Associated Eigenvectors." 
% ASME. J. Vib. Acoust. July 1996; 118(3): 390–397.
% https://doi.org/10.1115/1.2888195

clear; restoredefaultpath
%%

n_inner = 2;
n_ev    = 1;
degenerate = false;

x0 = 1/4; y0 = -1/4;

K_0 = eye(n_inner)*2;
dxK = ...
    [ 1 , 0 ;
      0 ,-1 ];
dyK = ...
    [ 0 , 1 ;
      1 , 0 ];

dx2K = zeros(n_inner);
dy2K = zeros(n_inner);
dx3K = zeros(n_inner);
dx4K = zeros(n_inner);
dy3K = zeros(n_inner);
dy4K = zeros(n_inner);

%%
demo