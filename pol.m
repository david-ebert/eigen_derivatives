function [Q0,polk] = pol(dL,k)
%POL Determine initial Polarization Q0 and decision matrix polk
%   [Q0,polk] = POL(dL)

n = length(dL{end});
d = length(dL);
if nargin == 1
    k = 1;
end
polk = ones(n)*k;

[Q0,dl] = eig(dL{k},"vector");

[dl,ind] = sort(dl);
Q0 = Q0(:,ind);
group = group_eig(dl);

Q_adj = eye(n);
for i = 1:max(group)
    I = group == i;
    if sum(I) > 1
        if k<d
            dL_sub = cell(d,1);
            for j = k+1:d; dL_sub{j} = Q0(:,I)'*dL{j}*Q0(:,I); end
            [Q_adj(I,I),polk(I,I)] = pol(dL_sub,k+1); % recursion
        else
            polk = Inf; % ran out of derivatives
        end
    end
end
Q0 = Q0*Q_adj;
if k==1 && any(polk(:)==Inf)
    warning('eig_deriv:pol_unseparable', ...
        ['No more derivatives of eigenvalues with respect to eigenspace available. ' ...
        'Stopped and assumed that eigenvalues are identical.'])
end
end