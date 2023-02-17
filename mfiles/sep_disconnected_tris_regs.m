%-------------------------------------------------------------------------------
%
% Separate disconnected triangle regions
%
%-------------------------------------------------------------------------------
function ts = sep_disconnected_tris_regs(p,t)

%---------------------------------------------------------------------------
disp('Get Matlab Triangulation Object')
TR     = triangulation(t, p);
disp('Get Local Neighbors')
neighs = neighbors(TR);
inans = find(isnan(neighs(:)) == 1);
iinfs = find(isinf(neighs(:)) == 1);
disp(['Num NaN: ',num2str(length(inans))])
disp(['Num Inf: ',num2str(length(iinfs))])

%-------------------------------------------------------------------------------
% Construct a connection matrix
disp('Make Sparse Adjancy Matrix')
% L = sparse(size(p,1),size(p,1));
% for i = 1:size(p,1)
%     %---------------------------------------------------------------------------
%     % Loop through all the neighboring nodes
%     for m = 1:3
%         j      = neighs(i,m);
%         L(i,j) = 1;
%         L(j,i) = 1;
%         L(i,i) = 1;
%     end
% end
spairs = zeros(size(t,1)*3*3,2);
k = 1;
for i = 1:size(t,1)
    %---------------------------------------------------------------------------
    % Loop through all the neighboring nodes
    for m = 1:3
        j      = neighs(i,m);
        if isnan(j) == 0
            spairs(k,:) = [i j];
            k = k+1;
            spairs(k,:) = [j i];
            k = k+1;
            spairs(k,:) = [i i];
            k = k+1;
        end
    end
end
spairs = spairs(1:k-1,:);
L = sparse(spairs(:,1),spairs(:,2),1+0*spairs(:,1),size(t,1),size(t,1),9*size(t,1));

%-------------------------------------------------------------------------------
disp('Find Local Connections')
Lcnx    = find_local_connections(L);
disp('Collect Connected Regions')
minLcnx = coll_cnxs(Lcnx);

%-------------------------------------------------------------------------------
for n = 1:length(minLcnx)
    ts{n} = t(minLcnx{n},:);
end

