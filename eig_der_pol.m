function [dL,du,Q0,polk,dQ,dl_Q] = eig_der_pol(l0,u0,A0,dA,M0,dM)
%EIG_DER_POL Calculate polarized derivatives
%   [dL,du,Q0,polk,dQ,dl_Q] = EIG_DER_POL(l0,u0,A0,dA) 
%   [dL,du,Q0,polk,dQ,dl_Q] = EIG_DER_POL(l0,u0,A0,dA,M0) 
%   [dL,du,Q0,polk,dQ,dl_Q] = EIG_DER_POL(l0,u0,A0,dA,M0,dM) 

switch nargin
    case 4
        [dL,du] = eig_der(l0(1),u0,A0,dA);
        [Q0,polk] = pol(dL);
        [dQ,dl_Q] = pol_der(dL,Q0,polk,u0,du);
    case 5
        [dL,du] = eig_der(l0(1),u0,A0,dA,M0);
        [Q0,polk] = pol(dL);
        [dQ,dl_Q] = pol_der(dL,Q0,polk,u0,du,M0);
    case 6
        [dL,du] = eig_der(l0(1),u0,A0,dA,M0,dM);
        [Q0,polk] = pol(dL);
        [dQ,dl_Q] = pol_der(dL,Q0,polk,u0,du,M0,dM);
end

end