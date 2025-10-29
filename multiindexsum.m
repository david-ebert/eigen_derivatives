function [list] = multiindexsum(o,n)
%MULTINDEXSUM all multiindces of specific length and order 
%   [list] = MULTIINDEXSUM(o,n)
%
% Input Arguments
%   o - Order of the multiindex
%       positive scalar
%   n - Length of the multiindex
%       positive scalar
%
% Output Arguments
%   list - List of multiindices
%       [nchoosek(n+o-1,n-1)]xn-matrix
%
%   See also MULTINOM

[C{1:n}] = ndgrid(0:o);
C = cell2mat(cellfun(@(x)(reshape(x,[],1)),C,'UniformOutput',false));
list = C(sum(C,2) == o,:);

end