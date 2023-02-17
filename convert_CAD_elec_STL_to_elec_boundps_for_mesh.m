%--------------------------------------------------------------------------
%
% This code takes an STL that describes a flat electrode array, extracts
% each electrodes boundary, and outputs these boundary points for other
% meshing codes. User input is needed for 1) the STL file name, 2) the
% orientation of the STL (perhaps?), and 3) the desired labeling of the
% electrodes. 
%
%--------------------------------------------------------------------------
clear
clc
close all

%--------------------------------------------------------------------------
% Add paths
addpath mfiles
addpath S:\digihisto\Ethan\software\Matlab_toolboxes\stlread

%--------------------------------------------------------------------------
% An STL file to process
[file,path] = uigetfile('*.STL');


%--------------------------------------------------------------------------
% load the STL file,
[vs,fs] = stlread([path,'\',file],1);


%--------------------------------------------------------------------------
% Determine the orientation of the STL - we want x & y to be in the
% electrode direction and z going through the electrodes
[nvec,cent,long_ax,shrt_ax] = get_nrmal_vec(vs);
nvec    = round(nvec);
long_ax = round(long_ax);
shrt_ax = round(shrt_ax);
nvec    = nvec/norm(nvec);
long_ax = long_ax/norm(long_ax);
shrt_ax = shrt_ax/norm(shrt_ax);
R       = [long_ax shrt_ax nvec];
vs      = (R*(vs'))';

%--------------------------------------------------------------------------
figure
trisurf(fs,vs(:,1),vs(:,2),vs(:,3))
lbl_fmt_fig('X (mm)','Y (mm)','','','Z (mm)',12)
axis equal

%--------------------------------------------------------------------------
% Get the top layer
tcs     = get_tcs(vs,fs);
is      = find( abs(tcs(:,3) - max(tcs(:,3))) < 1e-6);
fs      = fs(is,:);
[fs,vs] = remove_unused_nodes(fs,vs);
[fs,vs] = remove_repeated_nodes(fs,vs);
figure
trisurf(fs,vs(:,1),vs(:,2),vs(:,3))
view(2)

%--------------------------------------------------------------------------
% Separate the triangulations
ts  = sep_disconnected_tris_regs(vs,fs);
nel = length(ts);
%--------------------------------------------------------------------------
% For each separate triangle get the boundary. Then find the closest
% distance between any two boundary nodes and add other nodes to get a
% consistent spacing of points
figure
trisurf(fs,vs(:,1),vs(:,2),0*vs(:,3))
view(2)
hold on
for n = 1:nel
    TR   = triangulation(ts{n},vs(:,1:2));
    is   = unique(ts{n}(:));
    cent = mean(vs(is,1:2),1);
    vtmp = vs(is,1:2) - repmat(cent,length(is),1);
    thts = atan2(vtmp(:,2),vtmp(:,1));
    [tmp,sid] = sort(thts);
    iebds{n}  = is(sid); 
    
    %----------------------------------------------------------------------
    % Ge the minimum distance between points
    D = sqrt(comp_pairwise_distmat(vs(iebds{n},1:2)));
    D = D + diag(max(D(:))*ones(size(D,1),1));
    mind = min(D(:));
        
    %----------------------------------------------------------------------
    % Get a set of boundary points for each electrode
    nps = [];
    for k = 1:length(iebds{n})
        if k < length(iebds{n})
            k2 = k+1;
        else
            k2 = 1;
        end
        dst  = norm(vs(iebds{n}(k),:)-vs(iebds{n}(k2),:));
        nadd = round(dst/mind)+1;
        if nadd > 2
            ss = linspace(0,1,nadd)';ss = ss(1:end-1);
            nps = [nps; ...
                (1-ss)*vs(iebds{n}(k),1)+ss*vs(iebds{n}(k2),1) ...
                (1-ss)*vs(iebds{n}(k),2)+ss*vs(iebds{n}(k2),2)];
        else
            nps = [nps; vs(iebds{n}(k),1:2)];
        end
        
    end           
    %----------------------------------------------------------------------
    % If there is more than 100 points, then take every other.
    if size(nps,1) > 100
        nps = nps(1:2:end,:);
    end
    
    plot(nps(:,1),nps(:,2),'.','markersize',8)
    %----------------------------------------------------------------------
    ebnds{n} = nps;
    %----------------------------------------------------------------------
    % Ge the minimum distance between points
    D = sqrt(comp_pairwise_distmat(nps));
    D = D + diag(max(D(:))*ones(size(D,1),1));    
    eb_minds(n) = min(D(:));
    text(mean(nps(:,1)),mean(nps(:,2)),num2str(n),'color','red','fontsize',16)
    
end

%--------------------------------------------------------------------------
ecents = zeros(nel,2);
for n = 1:nel
    nps = ebnds{n};
    ecents(n,:) = mean(nps,1);
end

%--------------------------------------------------------------------------
% Reorder the electrodes
figure;set_fig_relsiz(0.6);
trisurf(fs,vs(:,1),vs(:,2),0*vs(:,3),'linestyle','none','facealpha',0.5)
view(2)
hold on
for n = 1:length(ebnds)
    nps = ebnds{n};
    plot(nps(:,1),nps(:,2),'.','markersize',8)
    % text(ecents(n,1),ecents(n,2),num2str(n),'color','red','fontsize',16)
end
axis equal;axis off
emap = zeros(nel,1);
for n = 1:nel
    % Click on the kth electrode
    title(['Click on the ',num2str(n),'th Electrode'],'fontsize',12)
    [x,y]   = ginput(1);
    [tmp,i] = min(sqrt(sum( (ecents - repmat([x y],size(ecents,1),1)).^2,2)  ));
    emap(n) = i;
    text(ecents(i,1),ecents(i,2),num2str(n))
    plot(x,y,'m')
    drawnow
end


%---------------------------------------------------------------------------
% Reorder the electrodes
ebnds  = ebnds(emap);
ecents = ecents(emap,:);
figure;set_fig_relsiz(0.6);
trisurf(fs,vs(:,1),vs(:,2),0*vs(:,3),'linestyle','none','facealpha',0.5)
view(2)
hold on
for n = 1:length(ebnds)
    nps = ebnds{n};
    plot(nps(:,1),nps(:,2),'.','markersize',8)
    text(ecents(n,1),ecents(n,2),num2str(n),'color','red','fontsize',16)
end
title('Double-check Order is correct in the saved data','fontsize',12)
axis equal;axis off

%---------------------------------------------------------------------------
eval(['save CAD_elecbounds/elec_bnds_',ifblank(file(1:end-4)),' ebnds'])




