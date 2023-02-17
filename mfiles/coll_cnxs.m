%-------------------------------------------------------------------------------
% Collect the connections
function [minLcnx,nrs] = coll_cnxs(Lcnx)

%-------------------------------------------------------------------------------
check           = 1;
n               = 1;
minLcnx{1}      = Lcnx{n};
f_ids           = minLcnx{n};
searched_ids    = zeros(size(Lcnx,1),1);
searched_ids(1) = n;
newf_ids        = [];
nrs             = 1;

while check == 1    
    %---------------------------------------------------------------------------
    % Loop through all the find index and collect the connected nodes
    oldsearch = minLcnx{nrs};
    for k = f_ids
        %-----------------------------------------------------------------------
        % Add not previously found connected nodes to search, next time
        newf_ids   = [newf_ids Lcnx{k}];
        %-----------------------------------------------------------------------
        % Add the new connected nodes
        minLcnx{nrs} = [minLcnx{nrs} Lcnx{k}];
        %-----------------------------------------------------------------------
        % Remove these nodes from the not yet search 
        searched_ids(k) = 1;
        
    end 
    %---------------------------------------------------------------------------
    % Reduce the new search indices to the unique ones that were not
    % previously searched
    newf_ids = setdiff(unique(newf_ids),oldsearch);
    
    %---------------------------------------------------------------------------
    % Set the search indices
    f_ids    = newf_ids;
    newf_ids = [];
    
    %---------------------------------------------------------------------------
    if isempty(f_ids) == 1    
        %-----------------------------------------------------------------------
        % Find the next unconnected node to start a new search
        newsearch_ids = find(searched_ids == 0);
        n             = min(newsearch_ids);
        
        %-----------------------------------------------------------------------
        % Display a connected region has been searched.
        nrs = nrs + 1;
        %         disp([num2str(nrs-1),' Regions have been searched'])
        %         disp([num2str(length(newsearch_ids)),' Remaining nodes'])
        %         disp(' ')
                

        %-----------------------------------------------------------------------
        % Check if we're done.
        if isempty(newsearch_ids) == 1
            %-------------------------------------------------------------------
            % If there are no more nodes then stop.
            check = 0;
        else
            %-------------------------------------------------------------------
            % Initialize the start (initial find indices and mark the node
            % searched)
            minLcnx{nrs}    = Lcnx{n};
            f_ids           = minLcnx{nrs};
            searched_ids(n) = 1;
        end
    end
end