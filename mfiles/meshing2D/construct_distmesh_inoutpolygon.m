%-------------------------------------------------------------------------------
%
% Construct a 2D mesh over an inner and outer polygon domain. The h-values
% of the mesh should match the edge spacing of the inner and outer
% boundaries. 
%
%-------------------------------------------------------------------------------
function [p,t] = construct_distmesh_inoutpolygon(pin,pout,dbg_flg,inchval,reorder)

if nargin < 5
    reorder = 1;
end
%-------------------------------------------------------------------------------
% Make sure the inner and outer boundaries are arranged in an CCW fashion
% Inner:
if reorder == 1
    ptmp     = pin - repmat(mean(pin,1),size(pin,1),1);
    ts       = atan2(ptmp(:,2),ptmp(:,1));
    [tmp,is] = sort(ts);
    pin      = pin(is,:);
    % Outer:
    ptmp     = pout - repmat(mean(pout,1),size(pout,1),1);
    ts       = atan2(ptmp(:,2),ptmp(:,1));
    [tmp,is] = sort(ts);
    pout     = pout(is,:);
end

%-------------------------------------------------------------------------------
% H-values: Calculate the spacing of points on the inside portion and the 
% outside portion
% Inner:
D      = comp_pairwise_distmat(pin);   
minDs  = sort(D,2,'ascend');
% hin    = 1.1 * sqrt(max( minDs(:,2))) * sin(60*pi/180);
if nargin < 4
    hin   = 1.25 * sqrt(max(minDs(:,2))) * sin(60*pi/180);
else
    hin   = inchval*1.25 * sqrt(max(minDs(:,2))) * sin(60*pi/180);
end
% Outer:
D      = comp_pairwise_distmat(pout);       
minDs  = sort(D,2,'ascend');
% hout   = 1.1 * sqrt(max(minDs(:,2))) * sin(60*pi/180);
if nargin < 4
    hout   = 1.25 * sqrt(max(minDs(:,2))) * sin(60*pi/180);
else
    hout   = inchval*1.25 * sqrt(max(minDs(:,2))) * sin(60*pi/180);
end
% Sort the two h-values (small to big)
hmin   = min([hin hout]);

%-------------------------------------------------------------------------------
% Distmesh: Make the 2D mesh with whole
%-------------------------------------------------------------------------------
% Construct a distance function with the drill bit cut out. In this
% function the drillbit is much more complicated than a circle
fd=@(p,binf) ddiff( dpoly(p,[pout; pout(1,:)]),dpoly(p,[pin; pin(1,:)]) );
%-------------------------------------------------------------------------------
% Set the bounding box and the fixed points
fpts = [pin; pout];
bbx  = [min(pout,[],1); max(pout,[],1)];    
%-------------------------------------------------------------------------------
% Construct the h-scaling function:
%   The scaling function is basic. For each input point we find the closest
%   Inner point and closest outer point. Then we construct a linear
%   function that changes from hin to hout based on where it is between the
%   two points. 
binf.hin  = hin;
binf.hout = hout;
binf.pin  = pin;
binf.pout = pout;

figure
[p,t]=distmesh2d(fd,@hscale_func,hmin,bbx,fpts,binf);%not relates to p above. p is n x 3, t is triangulation of dimen m x 3

%-------------------------------------------------------------------------------
if dbg_flg == 1 
    hold on
    plot(fpts(:,1),fpts(:,2),'dg','linewidth',2,'markersize',3)
    lbl_fmt_fig('X','Y','')
    axis equal
    
end

%-------------------------------------------------------------------------------
% Construct the h-scaling function:
%   The scaling function is basic. For each input point we find the closest
%   Inner point and closest outer point. Then we construct a linear
%   function that changes from hin to hout based on where it is between the
%   two points. 
function hsc = hscale_func(p,binf)

%-------------------------------------------------------------------------------
% Extract the boundary data
hin  = binf.hin;
hout = binf.hout;
pin  = binf.pin;
pout = binf.pout;

%-------------------------------------------------------------------------------
% Find the closest points for each p
closest = zeros(size(p,1),2);
for n = 1:size(p,1)
    [~,closest(n,1)] = min( ( pin(:,1)-p(n,1)).^2 + ( pin(:,2)-p(n,2)).^2 );
    [~,closest(n,2)] = min( (pout(:,1)-p(n,1)).^2 + (pout(:,2)-p(n,2)).^2 );
end

%-------------------------------------------------------------------------------
% Calculate the distance between the Inner to Outer closest points    
din2out = sqrt(sum( (pout(closest(:,2),:) - pin(closest(:,1),:)).^2,2));

%-------------------------------------------------------------------------------
% Construct a linear fit that scales between 1 and the ratio of the
% h-values
if hin < hout
    %---------------------------------------------------------------------------
    % Interior boundary spacing is smaller than Outer
    %---------------------------------------------------------------------------
    % Determine the distance of all ps to the closest interior points
    din2p   = sqrt(sum( ( p - pin(closest(:,1),:)).^2,2));
    
    %---------------------------------------------------------------------------
    % Calculate the scale for each: 
    % 
    %           scale = 1*(1-t) + hout/hin*t
    % 
    % where t is the distance from the inner boundary to p relative to the
    % outer boundary 
    t   = din2p ./ din2out;
    hsc = 1*(1-t) + hout/hin*t;    
        
else
    %---------------------------------------------------------------------------
    % Outer boundary spacing is smaller than Inner
    %---------------------------------------------------------------------------
    % Determine the distance of all ps to the closest outer points
    dout2p   = sqrt(sum( ( p - pin(closest(:,1),:)).^2,2));
    
    %---------------------------------------------------------------------------
    % Calculate the scale for each: 
    % 
    %           scale = 1*(1-t) + hout/hin*t
    % 
    % where t is the distance from the inner boundary to p relative to the
    % outer boundary 
    t   = dout2p ./ din2out;
    hsc = 1*(1-t) + hin/hout*t;
end

