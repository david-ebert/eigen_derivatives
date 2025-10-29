function [dL,du] = eig_der(l0,u0,A0,dA,M0,dM)
%EIG_DER - Derivatives of eigenpairs (with respect to the eigenspace)
%   [dL,du] = EIG_DER(l0,u0,A0,dA)
%   [dL,du] = EIG_DER(l0,u0,A0,dA,M0)
%   [dL,du] = EIG_DER(l0,u0,A0,dA,M0,dM)
% 
% Input Arguments
%   l0 - Eigenvalue of eigenspace
%       scalar
%   u0 - Eigenfunctions of eigenspace
%       nxm-matrix
%   A0 - (Stiffness) matrix
%       nxn-matrix
%   dA - Derivatives of (stiffness) matrix
%       cell vector of nxn-matrices
%   M0 - Mass matrix
%       eye(n) (default) nxn-matrix;
%   dM - Derivatives of mass matrix
%       (optional) cell vector of nxn-matrices
%       assumed zeros if not specified
%
% Output Arguments
%   dL - Eigenvalue derivatives wrt eigenspace
%       cell vector of mxm-matrices
%   du - Eigenfunction derivatives wrt eigenspace
%       cell vector of nxm-matrices

[n_dof,m] = size(u0);
d = length(dA);
L_dL = cell(d+1,1);
L_dL{1} = l0*eye(m);
u_du = cell(d+1,1);
u_du{1} = u0;
clear u0

A_dA = cell(d+1,1);
A_dA{1} = A0;
for k = 1:d
    A_dA{k+1} = dA{k};
end
clear A0 dA
switch nargin 
    case 4; M_dM{1} = speye(n_dof);
    case 5; M_dM{1} = M0;
    case 6
        M_dM = cell(d+1,1);
        M_dM{1} = M0;
        for k = 1:d
            M_dM{k+1} = dM{k};
        end
end
clear M0 dM

%% calc

nw_tile = A_dA{1} - M_dM{1}*l0;
ne_tile = -sparse(M_dM{1}*u_du{1});
lhs     = [nw_tile, ne_tile; ne_tile', sparse(m,m)];

for k = 1:d
    %% orthonormal conditions
    if nargin == 4
        diagterm = zeros(m);
    else
        diagterm  = zeros(m,1);
        mindexs = multiindexsum(k,3);
        mindexs = mindexs(mindexs(:,1)<k & mindexs(:,3)<k,:);
        if nargin < 6
            mindexs = mindexs(mindexs(:,2)==0,:);
        end
        mnc = multinom(mindexs);
        [n_build,~] = size(mindexs);
        for i_build = 1:n_build
            ind = mindexs(i_build,:);
            for i = 1:m
                diagterm(i) = diagterm(i) + ...
                    (mnc(i_build)/2)* u_du{1+ind(1)}(:,i)'*M_dM{1+ind(2)}*u_du{1+ind(3)}(:,i);
            end
        end
        diagterm = diag(diagterm);
    end
    %% main conditions
    rhs_n = zeros(n_dof,m);
    % A
    mindexs = multiindexsum(k,2);
    mindexs = mindexs(mindexs(:,2)<k,:);
    mnc = multinom(mindexs);
    [n_build,~] = size(mindexs);
    for i_build = 1:n_build
        ind = mindexs(i_build,:);
        rhs_n = rhs_n + ...
            - mnc(i_build)* A_dA{1+ind(1)} * u_du{1+ind(2)} ;
    end
    % M
    mindexs = multiindexsum(k,3);
    mindexs = mindexs(mindexs(:,2)<k & mindexs(:,3)<k,:);
    if nargin < 6
        mindexs = mindexs(mindexs(:,1)==0,:);
    end
    mnc = multinom(mindexs);
    [n_build,~] = size(mindexs);
    for i_build = 1:n_build
        ind = mindexs(i_build,:);
        rhs_n = rhs_n + ...
            + mnc(i_build)* M_dM{1+ind(1)} * u_du{1+ind(2)} * L_dL{1+ind(3)} ;
    end
    %% assemble & solve
    rhs = [rhs_n;diagterm];
    sol = lhs\rhs;
    u_du{1+k} = sol(1:n_dof,:);
    L_dL{1+k} = sol(n_dof+1:end,:);
end

%% output
dL = cell(d,1);
du = cell(d,1);
for k = 1:d
    dL{k} = L_dL{k+1};
    du{k} = u_du{k+1};
end

end
