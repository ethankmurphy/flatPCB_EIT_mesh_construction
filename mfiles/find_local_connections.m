%-------------------------------------------------------------------------------
% Find local connections in the regularization matrix
function Lcnx = find_local_connections(L,inc_self)

if nargin < 2
    inc_self = 0;    
end
%-------------------------------------------------------------------------------
% Initialize a cell array to contain all the connections
Lcnx = cell(size(L,1),1);

%-------------------------------------------------------------------------------
for n = 1:size(L,1)
    if inc_self == 0    
        Lcnx{n} = setdiff(find(abs(L(n,:)) > 0),n);
    else
        Lcnx{n} = find(abs(L(n,:)) > 0);
    end
end