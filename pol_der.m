function [dP,dl] = pol_der(dL,P0,polk,u0,du,M0,dM)
%POL_DER calculates available polarized derivatives of dl and corresponding dP
%   [dP,dl] = POL_DER(dL,P0,polk,u0,du)
%   [dP,dl] = POL_DER(dL,P0,polk,u0,du,M0)
%   [dP,dl] = POL_DER(dL,P0,polk,u0,du,M0,dM)
%
% Input Arguments
%   dL - Eigenvalue derivatives wrt eigenspace
%       cell vector of mxm-matrices
%   P0 - Initial polarization matrix
%       mxm-matrix
%   polk - Decision matrix
%       mxm-matrix
%   u0 - Eigenvectors of eigenspace
%       nxm-matrix
%   du - Eigenvector derivatives wrt eigenspace
%       cell vector of nxm-matrices
%   M0 - Mass matrix
%       eye(n) (default) nxn-matrix;
%   dM - Derivatives of mass matrix
%       (optional) cell vector of nxn-matrices
%       assumed zeros if not specified
%
% Output Arguments
%   dP - Derivatives of polarization matrix
%       cell vector of mxm-matrices
%   dl - Polarized eigenvalue derivatives
%       cell vector of mx1-vectors
%
% See also EIG_DER, POL, DU_POL

m = length(P0);
d = length(dL);

% merge u0 and du
u_du = cell(length(du)+1,1);
for i = 1:length(du)
    u_du{1+i} = du{i};
end
u_du{1} = u0;
clear u0 du

% merge M0 and dM
switch nargin 
    case 5; M_dM{1} = speye(size(u_du{1},1));
    case 6; M_dM{1} = M0;
    case 7
        M_dM = cell(d+1,1);
        M_dM{1} = M0;
        for k = 1:d
            M_dM{k+1} = dM{k};
        end
end
clear M0 dM

% merge P0 and initialize
P_dP = cell(d+1,1);
P_dP{1} = P0;
for i = 1:d
    P_dP{1+i} = zeros(m);
end
determined = false(d,m);
dl = cell(d,1);
for i = 1:d
    dl{i} = zeros(m,1);
end

%% calculations

for k_l = 1:d
    for i = 1:m
        if k_l <= polk(i,i) 
            dl{k_l}(i) = P0(:,i)'*dL{k_l}*P0(:,i);
            continue
        end
        lhs = -P0';
        rhs = zeros(m,1);
        kii = polk(i,i); % deciding order
        k_P = k_l-kii;

        %% adjust dP
        % set conditions related to previously determined dP 
        % this is necessary for calculation of dl
        for j = 1:m
            kij = polk(i,j);
            if kij < kii
                rhs(j) = coeff_builder(i,j,k_P,kij,P_dP,dL,dl);
            end
        end
        P_dP{1+k_P}(:,i) = lhs\rhs;

        %% dl
        % calculate dl separately
        sumR = dL{k_l}*P0(:,i);
        mindexs = multiindexsum(k_l,2);
        mindexs = mindexs(mindexs(:,1)>=kii & mindexs(:,1)<k_l,:); % dkL, k > polk
        mnc = multinom(mindexs);
        [n_build,~] = size(mindexs);
        for i_build = 1:n_build
            ind = mindexs(i_build,:);
            sumR = sumR + ...
                mnc(i_build) * (dL{ind(1)}*P_dP{1+ind(2)}(:,i) - P_dP{1+ind(2)}(:,i)*dl{ind(1)}(i));
        end
        dl{k_l}(i) = P0(:,i)'*sumR;

        %% dP final
        % use dl for condition and calculate dP
        for j = 1:m %% assemble condition to P0(:,j)
            kij = polk(i,j);
            if kij == kii
                if i~=j % standard conditions
                    rhs(j) = coeff_builder(i,j,k_P,kij,P_dP,dL,dl);
                else % normalization condition
                    mindexs = multiindexsum(k_P,5);
                    mindexs = mindexs(mindexs(:,1)<k_P & mindexs(:,5)<k_P,:);
                    if nargin < 7
                        mindexs = mindexs(mindexs(:,3)==0,:);
                    end
                    mnc_norm = multinom(mindexs);
                    [n_buildnorm,~] = size(mindexs);
                    for i_build = 1:n_buildnorm
                        ind = mindexs(i_build,:);
                        rhs(j) = rhs(j) + ...
                            (mnc_norm(i_build)/2) *...
                            P_dP{1+ind(1)}(:,i)'* (u_du{1+ind(2)}'*M_dM{1+ind(3)}*u_du{1+ind(4)}) *P_dP{1+ind(5)}(:,i);
                    end
                end
            end
        end
        P_dP{1+k_P}(:,i) = lhs\rhs;
        determined(k_P,i) = true;
    end  
end

%% output

dP = cell(length(P_dP)-1,1);
for i = 1:length(P_dP)-1
    dP{i} = P_dP{i+1};
    dP{i}(:,~determined(i,:)) = NaN;
end

end

function r = coeff_builder(i,j,k_P,kij,P_dP,dL,dl_P)
r = 0;
m = size(P_dP{1},1);
mindexs = multiindexsum(k_P + kij,2);
mindexs = mindexs(mindexs(:,1)>=kij,:); % dkL, k > polk
mnc = multinom(mindexs);
[n_build,~] = size(mindexs);
for i_build = 1:n_build
    ind = mindexs(i_build,:);
    r = r + ...
        mnc(i_build) * P_dP{1}(:,j)'* (dL{ind(1)}-dl_P{ind(1)}(i)*eye(m)) *P_dP{1+ind(2)}(:,i);
end
r = r/ (multinom([kij,k_P])* (dl_P{kij}(j) - dl_P{kij}(i)));
end