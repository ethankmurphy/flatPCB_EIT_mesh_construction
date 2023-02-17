%--------------------------------------------------------------------------
% 
% Construct the mesh using gmsh. The steps are the following:
%   1. Pick the electrode boundary to use
%   2. Construct the outersurface,
%   3. Run gmsh
%   4. Format the output of gmsh in our format, define the electrodes, and 
%      save the data
%       
% Other things are interest
%   * Defining the 'fine' region
%   * Defining the main mesh dimension
%   * Defining the h-values
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
%   2. distmesh - This is a matlab simple (mainly) 2D meshing software. 
%                 You just need to download it, https://github.com/ionhandshaker/distmesh
%                 and then add the path to it. 
addpath S:\digihisto\Ethan\software\distmesh

%--------------------------------------------------------------------------
% Pick the electrode boundaries of interest
[file,path] = uigetfile('*.MAT');

%--------------------------------------------------------------------------
% Construct the outersurface
fname = ['outersurface_',file(1:end-4),'_fac',ifdec(num2str(fac))];    
eval(['load dat/',fname,' p t elpts imain isd'])

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
% Run and save a mesh in gmsh
geo_str  = ['mesh_',file(1:end-4),'_h',ifdec(num2str(helec)), ...
    '_',ifdec(num2str(hfine)),'_',ifdec(num2str(hfar))];
eval(['load gmshes/',geo_str,' msh elpts'])



% Calculate the electrode areas
nds = msh.node;
for k = 1:length(msh.elec)
    ftmp      = msh.face(msh.elec{k},:);
    eareas(k) = sum(calc_TRI_area(ftmp,nds));
end
figure;hold on
plot(eareas*1000^2,'-ok')
lbl_fmt_fig('Electrode','Areas (mm^2)','','','',12)    

%---------------------------------------------
figure
plot_msh_elecs_only(msh,1,1,1)
% camlight left
axis equal
axis off
saveas(gcf,['figs/',geo_str,'_2D_downward_elecs_only'],'png')

%---------------------------------------------
figure;set(gcf,'position',[917         142        1696        1078])
hold on
simpplot(msh.node,msh.elem,'p(:,2)>0');hold on
trisurf(msh.face,msh.node(:,1),msh.node(:,2),msh.node(:,3),'linestyle','none','facecolor','blue','facealpha',0.3)
plot_msh_elecs_only(msh,1,0,1)
%---------------------------------------------
axis equal;axis off
view(3)
camlight left
saveas(gcf,['figs/',geo_str,'_3D_view'],'png')


msh


