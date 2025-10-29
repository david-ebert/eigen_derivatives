function [N] = multinom(alpha)
%MULTINOM multinomial coefficient of list of multiindeces
%   [N] = multinom(alpha)
%
% Input Arguments
%   alpha - List of multiindicies
%       nxm matrix
%
% Output Arguments
%   list - Multinomial coefficient of multiindices
%       nx1 vector
%
%   See also MULTIINDEXSUM

[n,l] = size(alpha); % length of multiindex

N = zeros(n,1);
for j = 1:n
    ord = sum(alpha(j,:)); % order of multiindex

    h = [];
    for i = 1:l
        h = [h , 1:alpha(j,i)]; %#ok<AGROW>
    end
    N(j) = prod((1:ord)./h);
end

end

