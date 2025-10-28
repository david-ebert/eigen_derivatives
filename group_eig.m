function [group] = group_eig(dl_Q,tol)
%GROUP_EIG Determines groups of eigenvalues
%   [group] = GROUP_EIG(dl_Q)
%   [group] = GROUP_EIG(dl_Q,tol)

if nargin == 1
    tol = 1E-5;
end

n = length(dl_Q);
group = zeros(n,1);
group(1) = 1;
startofgroup = 1;
numberofgroup = 1;
for i=2:n
    if (dl_Q(i) - dl_Q(startofgroup)) < tol
        group(i) = numberofgroup;
    else
        startofgroup = i;
        numberofgroup = numberofgroup +1;
        group(i) = numberofgroup;
    end
end

end