function [du_Q,u0_Q] = du_pol(u0,du,Q0,dQ)
%DU_POL polarized eigenfunction derivatives
%   [du_Q] = DU_POL(u0,du,Q0):     constant polarization
%   [du_Q] = DU_POL(u0,du,Q0,dQ):  non-constant polarization

u0_Q = u0*Q0;
switch nargin 
    case 3
        d = length(du);
        du_Q = cells(d,1);
        for i = 1:d
            du_Q = du{1}*Q0;
        end
    case 4
        d = length(dQ);
        size_du = size(u0);
    
        % merge u0 and du
        u_du = cell(length(du)+1,1);
        for i = 1:length(du)
            u_du{1+i} = du{i};
        end
        u_du{1} = u0;
        clear u0 du
    
        % merge Q0 and dQ
        Q_dQ = cell(d+1,1);
        Q_dQ{1} = Q0;
        for i = 1:d
            Q_dQ{1+i} = dQ{i};
        end
        clear Q0 dQ
    
        %% calculations
        du_Q = cell(d,1);
        for k = 1:d
            mindexs = multiindexsum(k,2);
            mnc = multinom(mindexs);
            [n_build,~] = size(mindexs);
    
            dku = zeros(size_du);
            for i_build = 1:n_build
                ind = mindexs(i_build,:);
                dku = mnc(i_build) * u_du{1+ind(1)}* Q_dQ{1+ind(2)};
            end
            du_Q{k} = dku;
        end
end

end
