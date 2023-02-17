%-------------------------------------------------------------------------------
%
% Get the center of each electrode. This assumes a standard NDRM mesh,
% which is defined as a structure array of the following form
% 
% msh.node - Nodes of a FEM mesh, size N_n x 3
% msh.elem - Tetrahedra elements of FEM mesh, size N_e x 4
% msh.face - Face of elements of the FEM mesh (generally only boundary
%            faces), sice N_f x 3.
% msh.elec - A N_el x 1 cell array where N_el is the number of electrodes.
%            Each cell is size 1 x N_elf_i where N_elf_i is the number of
%            face elements of the ith electrode. 
%
%-------------------------------------------------------------------------------
function ecents = get_el_cents(msh)

Nel    = length(msh.elec);
ecents = zeros(Nel,3);
fcs    = get_tcs(msh);
for k = 1:Nel
    ecents(k,:) = mean(fcs(msh.elec{k},:),1);    
    % ecents(k,:) = min(fcs(msh.elec{k},:),[],1);    
end