%-------------------------------------------------------------------------------
%
% Write the boundary points to a geo file
%
%------------------------------------------------------------------------------- 
function write_compoundvol2geo(prfx,vol_inds1,vol_inds2,ios)   

%-------------------------------------------------------------------------------
% Open the geo file
fid = fopen([prfx,'_vol.geo'],'w');

%-------------------------------------------------------------------------------
% Construct 2 volumes from the set of triangle surfaces
%-------------------------------------------------------------------------------
% First one
fprintf(fid,'Surface Loop(%i) = {',1+ios);
for n = 1:length(vol_inds1)-1
    fprintf(fid,'%i,',vol_inds1(n)+ios);
end
fprintf(fid,'%i};\n',vol_inds1(end)+ios);
%-------------------------------------------------------------------------------
fprintf(fid, ...
    'Volume(%i) = {%i};\n'              , ...
    1+ios                               , ...
    1+ios );
%-------------------------------------------------------------------------------
% Second one
fprintf(fid,'Surface Loop(%i) = {',2+ios);
for n = 1:length(vol_inds2)-1
    fprintf(fid,'%i,',vol_inds2(n)+ios);
end
fprintf(fid,'%i};\n',vol_inds2(end)+ios);
%-------------------------------------------------------------------------------
fprintf(fid, ...
    'Volume(%i) = {%i};\n'              , ...
    2+ios                               , ...
    2+ios );

%-------------------------------------------------------------------------------
% Combine the volumes using the Compound Volume function
fprintf(fid, ...
    'Compound Volume (%i) = {%i, %i};\n', ...
    3+ios                               , ...
    1+ios                               , ...
    2+ios);

%-------------------------------------------------------------------------------
fclose(fid);