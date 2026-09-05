function [du_pol,u0_pol] = du_pol(u0,du,P0,dP)
%DU_POL - Polarized eigenvector derivatives
%   [du_P] = DU_POL(u0,du,P0):     constant polarization
%   [du_P] = DU_POL(u0,du,P0,dP):  non-constant polarization
%
% Input Arguments
%   u0 - Initial eigenvectors
%       nxm-matrix
%   du - Derivatives of eigenvectors wrt eigenspace
%       cell vector of nxm-matrices
%   P0 - Initial polarization matrix
%       mxm-matrix
%   dP - Derivatives of polarization matrix
%       cell vector of mxm-matrices
%
% Output Arguments
%   du_pol - Polarized eigenvector derivatives
%       cell vector of nxm-matrices
%   u0_pol - Polarized initial eigenvectors
%       nxm-matrix
%
% See also EIG_DER, POL, POL_DER

u0_pol = u0*P0;
switch nargin 
    case 3
        d = length(du);
        du_pol = cell(d,1);
        for i = 1:d
            du_pol{i} = du{i}*P0;
        end
    case 4
        d = length(dP);
        size_du = size(u0);
    
        % merge u0 and du
        u_du = cell(length(du)+1,1);
        for i = 1:length(du)
            u_du{1+i} = du{i};
        end
        u_du{1} = u0;
        clear u0 du
    
        % merge P0 and dP
        P_dP = cell(d+1,1);
        P_dP{1} = P0;
        for i = 1:d
            P_dP{1+i} = dP{i};
        end
        clear P0 dP
    
        %% calculations
        du_pol = cell(d,1);
        for k = 1:d
            mindexs = multiindexsum(k,2);
            mnc = multinom(mindexs);
            [n_build,~] = size(mindexs);
    
            dku = zeros(size_du);
            for i_build = 1:n_build
                ind = mindexs(i_build,:);
                dku = dku + mnc(i_build) * u_du{1+ind(1)}* P_dP{1+ind(2)};
            end
            du_pol{k} = dku;
        end
end

end
