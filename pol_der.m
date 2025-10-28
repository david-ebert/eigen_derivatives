function [dQ,dl_Q] = pol_der(dL,Q0,polk,u0,du,M0,dM)
%POL_DER calculates available polarized derivatives of dl and corresponding dQ
%   [dQ,dl_Q] = POL_DER(dL,Q0,polk,u0,du)
%   [dQ,dl_Q] = POL_DER(dL,Q0,polk,u0,du,M0)
%   [dQ,dl_Q] = POL_DER(dL,Q0,polk,u0,du,M0,dM)

m = length(Q0);
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

% merge Q0 and initalize
Q_dQ = cell(d+1,1);
Q_dQ{1} = Q0;
for i = 1:d
    Q_dQ{1+i} = zeros(m);
end
dl_Q = cell(d,1);
for i = 1:d
    dl_Q{i} = zeros(m,1);
end

%% calculations

for k_l = 1:d
    for i = 1:m
        if k_l <= polk(i,i) 
            dl_Q{k_l}(i) = Q0(:,i)'*dL{k_l}*Q0(:,i);
            continue
        end
        lhs = -Q0';
        rhs = zeros(m,1);
        kii = polk(i,i); % deciding order
        k_Q = k_l-kii;

        %% adjust dQ
        % set conditions related to previously determined dQ 
        % this is necessary for calculation of dl
        for j = 1:m
            kij = polk(i,j);
            if kij < kii
                rhs(j) = coeff_builder(i,j,k_Q,kij,Q_dQ,dL,dl_Q);
            end
        end
        Q_dQ{1+k_Q}(:,i) = lhs\rhs;

        %% dl
        % calculate dl separately
        sumR = dL{k_l}*Q0(:,i);
        mindexs = multiindexsum(k_l,2);
        mindexs = mindexs(mindexs(:,1)>=kii & mindexs(:,1)<k_l,:); % dkL, k > polk
        mnc = multinom(mindexs);
        [n_build,~] = size(mindexs);
        for i_build = 1:n_build
            ind = mindexs(i_build,:);
            sumR = sumR + ...
                mnc(i_build) * (dL{ind(1)}*Q_dQ{1+ind(2)}(:,i) - Q_dQ{1+ind(2)}(:,i)*dl_Q{ind(1)}(i));
        end
        dl_Q{k_l}(i) = Q0(:,i)'*sumR;

        %% dQ final
        % use dl for condition and calculate dQ
        for j = 1:m %% assemble condition to Q0(:,j)
            kij = polk(i,j);
            if kij == kii
                if i~=j % standard conditions
                    rhs(j) = coeff_builder(i,j,k_Q,kij,Q_dQ,dL,dl_Q);
                else % normalization condition
                    mindexs = multiindexsum(k_Q,5);
                    mindexs = mindexs(mindexs(:,1)<k_Q & mindexs(:,5)<k_Q,:);
                    if nargin < 7
                        mindexs = mindexs(mindexs(:,3)==0,:);
                    end
                    mnc_norm = multinom(mindexs);
                    [n_buildnorm,~] = size(mindexs);
                    for i_build = 1:n_buildnorm
                        ind = mindexs(i_build,:);
                        rhs(j) = rhs(j) + ...
                            (mnc_norm(i_build)/2) *...
                            Q_dQ{1+ind(1)}(:,i)'* (u_du{1+ind(2)}'*M_dM{1+ind(3)}*u_du{1+ind(4)}) *Q_dQ{1+ind(5)}(:,i);
                    end
                end
            end
        end
        Q_dQ{1+k_Q}(:,i) = lhs\rhs;
    end  
end

%% output

dQ = cell(length(Q_dQ)-1,1);
for i = 1:length(Q_dQ)-1
    dQ{i} = Q_dQ{i+1};
end

end

function r = coeff_builder(i,j,k_Q,kij,Q_dQ,dL,dl_Q)
r = 0;
m = size(Q_dQ{1},1);
mindexs = multiindexsum(k_Q + kij,2);
mindexs = mindexs(mindexs(:,1)>=kij,:); % dkL, k > polk
mnc = multinom(mindexs);
[n_build,~] = size(mindexs);
for i_build = 1:n_build
    ind = mindexs(i_build,:);
    r = r + ...
        mnc(i_build) * Q_dQ{1}(:,j)'* (dL{ind(1)}-dl_Q{ind(1)}(i)*eye(m)) *Q_dQ{1+ind(2)}(:,i);
end
r = r/ (multinom([kij,k_Q])* (dl_Q{kij}(j) - dl_Q{kij}(i)));
end