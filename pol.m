function [P0,polk] = pol(dL,tol,k)
%POL - Initial Polarization P0 and decision matrix polk
%   [P0,polk] = POL(dL)
%   [P0,polk] = POL(dL,tol)
%
% Input Arguments
%   dL - Eigenvalue derivatives wrt eigenspace
%       cell vector of mxm-matrices
%   tol - Tolerance for degeneracy
%       1E-5 (default) positive scalar
%
% Output Arguments
%   P0 - Initial polarization matrix
%       mxm-matrix
%   polk - Decision matrix
%       mxm-matrix
%
% See also EIG_DER, POL_DER

switch nargin 
    case 1
        k = 1;
        tol = 1E-5;
    case 2
        k = 1;
end
n = length(dL{end});
d = length(dL);

polk = ones(n)*k;

[P0,dl] = eig(dL{k},"vector");

[dl,ind] = sort(dl);
P0 = P0(:,ind);
group = group_eig(dl,tol);

P_adj = eye(n);
for i = 1:max(group)
    I = group == i;
    if sum(I) > 1
        if k<d
            dL_sub = cell(d,1);
            for j = k+1:d; dL_sub{j} = P0(:,I)'*dL{j}*P0(:,I); end
            [P_adj(I,I),polk(I,I)] = pol(dL_sub,tol,k+1); % recursion
        else
            polk = Inf(n); % ran out of derivatives
        end
    end
end
P0 = P0*P_adj;
if k==1 && any(polk(:)==Inf)
    warning('eig_deriv:pol_unseparable', ...
        ['No more derivatives of eigenvalues with respect to eigenspace available. ' ...
        'Stopped and assumed that eigenvalues are identical.'])
end
end