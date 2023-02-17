%-------------------------------------------------------------------------------
%
% Write the boundary points to a geo file
%
%------------------------------------------------------------------------------- 
function write_nodes2geo(prfx,pts,hvals,ios)   

%-------------------------------------------------------------------------------
% Open the geo file
fid = fopen([prfx,'_nodes.geo'],'w');
fprintf(fid,'p    = %i;\n',ios);
%-------------------------------------------------------------------------------
% Loop through the points
for n = 1:size(pts,1)
    fprintf(fid,'p = p+1;\n');
    fprintf(fid, ...
        'Point(p) = {%12.12f,%12.12f,%12.12f,%12.12f};\n', ...
        pts(n,1)        , ...
        pts(n,2)        , ...
        pts(n,3), ...
        hvals(n));
end


%-------------------------------------------------------------------------------
% Close the file
fclose(fid);