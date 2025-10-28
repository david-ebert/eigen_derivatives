function [list] = multiindexsum(order,n)
%MULTINDEXSUM all multiindces of lenght n with order 

[C{1:n}] = ndgrid(0:order);
C = cell2mat(cellfun(@(x)(reshape(x,[],1)),C,'UniformOutput',false));
list = C(sum(C,2) == order,:);

end
