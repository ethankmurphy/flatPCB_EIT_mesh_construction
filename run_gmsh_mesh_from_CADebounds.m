%--------------------------------------------------------------------------
% 
% Construct the mesh using gmsh. The steps are the following:
%   1. Pick the electrode boundary to use
%   2. Construct the outersurface,
%   3. Set the h-values for the nodes in the surface mesh: on/near
%      electrodes, in the subdomain, and in the outer region
%   3. Run gmsh
%   4. Format the output of gmsh in our format, define the electrodes, and 
%      save the data
%       
% There are some heuristic parameters that are automatically chosen. This 
% could be adjusted if needed. These parameters control 
%   1) the shape of the subdomain, 
%   2) the size of the subdomain and full domain, and 
%   3) the h-values on the electrodes, subdomain, and outer region
% 
% Here are some details of the decisions:
%   1) The shape of the subdomain: The shape of the domain is given by a
%   best-fitting ellipse (covariance of the distribution of CAD electrode
%   boundary points)
%
%   2) The size of the subomdain: We define it as a factor times the
%   electrode array boundary. So, the subdomain follows the shape from
%   above while always being at least the fac x the largest extent of the
%   electrodes in any direction. In the z-direction we choose the largest
%   electrode extent times the factor down in the z-direction. 
%
%   3) The full domain size: The full domain is defined by +/-Nl,+/-Nl, and
%   [-Nl, 0], where Nl is 5x(the largest extent of the subdomain top
%   plane).
%
%   4) The electrode h-values are 0.5 times the minimum of the maximum 
%      distances between any electrode. For instance, if all the electrodes
%      are circles, then this says the h-value is the radius.
% 
%   5) The subdomain h-values are equal to the the minimum of the maximum 
%      distances between any electrode. In the above example, this means
%      the h-value is the diameter of the electrodes.
% 
%   6) The h-values for the outer region is 1/4 the maximum distance
%      between points defining the subdomain. 
% 
%--------------------------------------------------------------------------
clear
close all
clc

%--------------------------------------------------------------------------
fac      = 2; % Factor to increase the subdomain about the electrode array


%--------------------------------------------------------------------------
% Set the path
addpath(genpath('mfiles'))
%---------------------------------------------
% Required software
%   1. gmsh - below is a string to the executable, notice the space after
%             the gmsh in the string
% rungmsh_str = '/jumbo/digihisto/Ethan/software/gmsh-4.3.0-Linux64/bin/gmsh ';
rungmsh_str = 'C:\gmsh-4.8.4-Windows64\gmsh ';
%---------------------------------------------
%   2. distmesh - This is a matlab simple (mainly) 2D meshing software. 
%                 You just need to download it, https://github.com/ionhandshaker/distmesh
%                 and then add the path to it. 
addpath S:\digihisto\Ethan\software\distmesh

%--------------------------------------------------------------------------
% Pick the electrode boundaries of interest
% [file,path] = uigetfile('*.MAT');
path = 'S:\digihisto\Ethan\Githubs\mesh_flat_elec_array_from_STL\CAD_elecbounds\';
% path = '/jumbo/digihisto/Ethan/Githubs/mesh_flat_elec_array_from_STL/CAD_elecbounds/';
file = 'elec_bnds_Electrode_Array_8ch_9o5mm_v1.mat'
% file = 'elec_bnds_sarcopenia_US_EIT_Electrode_Array.mat'

%--------------------------------------------------------------------------
% Construct the outersurface
fname = ['outersurface_',file(1:end-4),'_fac',ifdec(num2str(fac))];
if filechecker('dat',[fname,'.mat']) == 0
    [p,t,elpts,imain,isd] = construct_outer_surf(path,file,fac,1);
    eval(['save dat/',fname,' p t elpts imain isd'])
else
    eval(['load dat/',fname,' p t elpts imain isd'])
end
%--------------------------------------------------------------------------
% Set the h-values, let's make it a function of the maximum distance of 
%   * electrode points (i.e. max distance between any points describing a
%     single electrode)
%   * the subdomain
%   * the entire domain
% Electrodes: Pick the min of the maximum spacing
maxels = zeros(length(elpts),1);
for n = 1:length(elpts)
    D         = sqrt(comp_pairwise_distmat(elpts{n}));
    maxels(n) = max(D(:));
end
% Subdomain
tcs   = get_tcs(p,t);
maxsd = max(max(sqrt(comp_pairwise_distmat(tcs(isd,:)))));


%--------------------------------------------------------------------------
% Set parameters
helec     = round(min(maxels)/4,1) %  0.2           % // h-value on electrodes
hfine     = round(min(maxels)/2,1) % 0.4           % // h-value in the fine area
hfar      = round(maxsd/4,1) % 5

%--------------------------------------------------------------------------
% Set the h-values: Nodes 
hvals = hfar*ones(size(p,1),1);
% First set the fine nodes to be nodes associated with the subdomain
for n = 1:3
    hvals(t(isd,n)) = hfine;
end
% Second for all nodes within 1 mm of every electrode set it to an
% electrode h-value
for n = 1:length(elpts)
    for k = 1:length(elpts{n})
        is = find( sqrt(sum( (p - repmat([elpts{n}(k,:) 0],size(p,1),1)).^2,2)) < maxels(n));
        hvals(is) = helec;
    end
end
i1s = find( hvals == hfar);
i2s = find( hvals == hfine);
i3s = find( hvals == helec);
figure
hold on
plot3(p(i1s,1),p(i1s,2),p(i1s,3),'.k','markersize',8)
plot3(p(i2s,1),p(i2s,2),p(i2s,3),'.r','markersize',8)
plot3(p(i3s,1),p(i3s,2),p(i3s,3),'.g','markersize',8)


%--------------------------------------------------------------------------
% Construct the gmsh mesh
construct_gmsh_compound_mesh(p,t,hvals,imain,isd,rungmsh_str);

%--------------------------------------------------------------------------
% 2. Extract the nodes and elements from the gmsh mesh file
msh = read_gmsh_mesh_v3('gmshes/tmp_mesh',0);

%--------------------------------------------------------------------------
% 2.b.  Determine the face elements using NDRM's face mex face
%     finding algorithm
disp('Finding faces')
TR       = triangulation(msh.elem, msh.node);
msh.face = freeBoundary(TR);

%--------------------------------------------------------------------------
% Define the electrodes
tcs  = get_tcs(msh.node,msh.face);
itop = find(abs(tcs(:,3) - 0)<1e-3);
for n = 1:length(elpts)
    % ks = convhull(elpts{n}(:,1),elpts{n}(:,2));
    In = inpolygon(tcs(itop,1),tcs(itop,2),elpts{n}(:,1),elpts{n}(:,2));
    elec{n,1} = itop(find( In ==1))';
end
msh.elec = elec;
msh.node = msh.node/1000;

%--------------------------------------------------------------------------
% Run and save a mesh in gmsh
geo_str  = ['mesh_',file(1:end-4),'_h',ifdec(num2str(helec)), ...
    '_',ifdec(num2str(hfine)),'_',ifdec(num2str(hfar))];
eval(['save gmshes/',geo_str,' msh elpts'])
