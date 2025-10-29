function [group] = group_eig(l,tol)
%GROUP_EIG Groups of eigenvalues according to degeneracy
%   [group] = GROUP_EIG(l)
%   [group] = GROUP_EIG(l,tol)
%
% Input Arguments
%   l - Eigenvalues 
%       sorted vector (min to max)
%   tol - Tolerance for degeneracy
%       1E-5 (default) positive scalar
%
% Output Arguments
%   group - Index of Group
%       vector of integers

if nargin == 1
    tol = 1E-5;
end

n = length(l);
group = zeros(n,1);
group(1) = 1;
startofgroup = 1;
numberofgroup = 1;
for i=2:n
    if (l(i) - l(startofgroup)) < tol
        group(i) = numberofgroup;
    else
        startofgroup = i;
        numberofgroup = numberofgroup +1;
        group(i) = numberofgroup;
    end
end

end