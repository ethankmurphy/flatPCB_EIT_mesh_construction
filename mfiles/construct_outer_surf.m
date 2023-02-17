%--------------------------------------------------------------------------
%
%
%
%--------------------------------------------------------------------------
function [p,t,elpts,imain,isd] = construct_outer_surf(path,fname,fac,dbg_flg)

%--------------------------------------------------------------------------
% Load the mat-file of the electrode boundaries
eval(['load ',path,fname])
fps = [];
for n = 1:length(ebnds)
    fps   = [fps; ebnds{n}];
end
cent  = mean(fps,1);
fps   = fps - repmat(cent,size(fps,1),1);
% Demean the electrode points (if it was not centered
if norm(cent) > 1e-6
    for n = 1:length(ebnds)
        ebnds{n} = ebnds{n} - repmat(cent,size(ebnds{n},1),1);
    end
end
elpts = ebnds;

% R0  = Rotz(2*pi*rand);
% fps = (R0(1:2,1:2)*(fps'))';

%--------------------------------------------------------------------------
% Get boundary of fine area
[a,b,t0] = get_bound_fine_ell(fps,fac,dbg_flg);


%--------------------------------------------------------------------------
% Construct the intermediate domain mesh
zbnd    = [-max(sqrt(sum(fps.^2,2))*fac) 0];
perim   = pi*(3*(a+b) - sqrt((3*a+b)*(a+3*b)));
nr      = 64;
zdel    = abs(zbnd(1))/5;
[cylpts,cyltri] = construct_ellcyl_nobottop(a,b,t0,zbnd,nr,zdel,1);
is              = find(abs(cylpts(:,3)-min(cylpts(:,3)))<1e-3);
obs             = cylpts(is,1:2);
[cylb_p,cylb_t] = construct_distmesh_polygon_nofixed(obs,0);
cylb_p          = [cylb_p zbnd(1)+0*cylb_p(:,1)];
[cylpts,cyltri] = combine_triangulations(cylb_p,cylpts,cylb_t,cyltri,0);

%--------------------------------------------------------------------------
% Get the outer boundary, defined as the top of the intermediate
% boundary
is  = find( abs( cylpts(:,3) - 0 ) < 1e-3);
obs = cylpts(is,1:2);
if dbg_flg == 1
    figure
    hold on
    plot(obs(:,1),obs(:,2),'.k','markersize',8)
    plot(fps(:,1),fps(:,2),'.r','markersize',8)
    axis equal
end
savfnam = ['topsrfmsh_',fname(1:end-4),'_fac',ifdec(num2str(fac))];
if filechecker('dat',[savfnam,'.mat']) == 0   
    [ptop,ttop] = construct_distmesh_polygon_innerfixed(obs,fps,0,1.2);
    ptop        = [ptop 0*ptop(:,1)];
    eval(['save dat/',savfnam,' ptop ttop elpts'])
else
    eval(['load dat/',savfnam,' ptop ttop elpts'])
end

if dbg_flg == 1
    figure
    trisurf(ttop,ptop(:,1),ptop(:,2),ptop(:,3),'facealpha',0.5)
    hold on
    plot(fps(:,1),fps(:,2),'.r','markersize',8)
    plot(obs(:,1),obs(:,2),'.m','markersize',8)
    axis equal
    lbl_fmt_fig('X (cm)','Y (cm)','2D Distmesh Top Mesh with Inclusion','','',12)
    view(2)
    % axis(4*[-1 1 -1 1])
    % saveas(gcf,'figs/distmesh_topelecs','png')
end

%--------------------------------------------------------------------------
% Construct the total top mesh outer boundary
maxdst   = max(sqrt(sum(ptop.^2,2)));
Nlen_far = 5*maxdst;
xs       = linspace(-Nlen_far,Nlen_far,10)';
obs      = [ ...
    xs -Nlen_far+0*xs; ...
    xs  Nlen_far+0*xs; ...
    Nlen_far+0*xs  xs; ...
    -Nlen_far+0*xs xs; ...
    ];
obs = unique(obs,'rows');
is  = find( abs( cylpts(:,3) - 0 ) < 1e-3);
ibs = cylpts(is,1:2);
whos
[ptmp,ttmp] = construct_distmesh_inoutpolygon(ibs,obs,dbg_flg,1.2);
ptmp        = [ptmp 0*ptmp(:,1)];
[ptop,ttop] = combine_triangulations(ptop,ptmp,ttop,ttmp,0);
[ptop,ttop] = combine_triangulations(ptop,cylpts,ttop,cyltri,0);

%--------------------------------------------------------------------------
% Make the other sides of the box
%--------------------------------------------------------------------------
% Construct the outer boundary
xs       = linspace(-Nlen_far,Nlen_far,10)';
zs       = linspace(-Nlen_far,0,5)';
obs      = [ ...
    xs -Nlen_far+0*xs; ...
    xs  0*xs; ...
    Nlen_far+0*zs  zs; ...
    -Nlen_far+0*zs zs; ...
    ];
obs = unique(obs,'rows');
[ptemp,ttemp] = construct_distmesh_polygon_nofixed(obs,0);

psid1 = [-Nlen_far+0*ptemp(:,1) ptemp(:,1)  ptemp(:,2)];
psid2 = [Nlen_far+0*ptemp(:,1)  ptemp(:,1)  ptemp(:,2)];
psid3 = [ptemp(:,1) -Nlen_far+0*ptemp(:,1)  ptemp(:,2)];
psid4 = [ptemp(:,1)  Nlen_far+0*ptemp(:,1)  ptemp(:,2)];
obs      = [ ...
    xs -Nlen_far+0*xs; ...
    xs  Nlen_far+0*xs; ...
    Nlen_far+0*xs  xs; ...
    -Nlen_far+0*xs xs; ...
    ];
obs = unique(obs,'rows');
[ptemp,ttemp2] = construct_distmesh_polygon_nofixed(obs,0);
pbot  = [ptemp(:,1)  ptemp(:,2)       -Nlen_far+0*ptemp(:,1)];
%-------------------------------------------------------------------------------
% Combine triangulations
[p,t] = combine_triangulations(ptop,psid1,ttop,ttemp,0);
[p,t] = combine_triangulations(   p,psid2,   t,ttemp,0);
[p,t] = combine_triangulations(   p,psid3,   t,ttemp,0);
[p,t] = combine_triangulations(   p,psid4,   t,ttemp,0);
[p,t] = combine_triangulations(   p,pbot,    t,ttemp2,1);

%-------------------------------------------------------------------------------
% Get the element indices for the main domain and the inner domain
% The inner domain includes the inner top and the cylinder, whereas the
% main domain includes everything except the inner top.
tcs  = get_tcs(p,t);
ts   = atan2(tcs(:,2),tcs(:,1));
rell = a*b*sqrt( 1 ./ (b.^2.*cos(ts-t0).^2 + a.^2.*sin(ts-t0).^2));
rtcs = sqrt(sum(tcs(:,1:2).^2,2));
% itop = find( (sqrt(sum(tcs(:,1:2).^2,2)) < rad_mid) & (abs(tcs(:,3)-0) < 1e-3) );
% iinr = find( (sqrt(sum(tcs(:,1:2).^2,2)) < 2*rad_mid) & (tcs(:,3) < -1e-3) & (tcs(:,3) > -5-1e-3) );
itop = find( (rtcs < rell) & (abs(tcs(:,3)-0) < 1e-3) );
iinr = find( (rtcs < 2*rell) & (tcs(:,3) < -1e-3) & (tcs(:,3) > zbnd(1)-1e-3) );
iotr = setdiff(setdiff((1:size(t,1))',itop),iinr);
imain = [iotr; iinr];
isd   = [itop; iinr];

if dbg_flg == 1

    figure
    subplot(1,2,1)
    trisurf(t(imain,:),p(:,1),p(:,2),p(:,3),'facecolor','cyan','facealpha',0.5)
    title('main domain')
    axis equal
    subplot(1,2,2)
    trisurf(t(isd,:),p(:,1),p(:,2),p(:,3),'facecolor','cyan','facealpha',0.5)
    title('subdomain')
    axis equal
end


