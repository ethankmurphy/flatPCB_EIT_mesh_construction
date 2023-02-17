%--------------------------------------------------------------------------
%
% Construct an elliptic cylindrical mesh with no bottom or top
% The eleptical cylinder is defined by a major, minor axis, and offset
% angle (direction of major axis), which are assumed to be in the
% x/y-direction. The z-direction is assumed to be the height of the
% cylinder. 
%
%--------------------------------------------------------------------------
function [cylpts,cyltri] = construct_ellcyl_nobottop(a,b,t0,zbnd,nr,zdel,plot_flg)

%--------------------------------------------------------------------------
if nargin == 2
    nr   = 32;
    zdel = 1;
elseif nargin == 3    
    zdel = 1;
end

%--------------------------------------------------------------------------
% Get the z values
nz     = max([ (zbnd(2)-zbnd(1))/zdel 3]);
zs     = linspace(zbnd(1),zbnd(2),nz );

%--------------------------------------------------------------------------
% Get ellipse angles that go around the ellipse in equal steps
perim_approx = pi*(3*(a+b) - sqrt((3*a+b)*(a+3*b)))
perim   = quad(@(theta) polar_integrand(theta,a,b,t0),0,2*pi);
ds      = perim/(nr);
ts      = zeros(nr,1);
for n = 2:nr
    bestt = fminbnd( @(theta) abs(ds*(n-1) - quad(@(theta) polar_integrand(theta,a,b,t0),0,theta) ),0.001,2*pi ); 
    ts(n) = bestt;
end

% ts      = linspace(0,2*pi,1000)';ts = ts(1:end-1);
% rell    = a*b*sqrt( 1 ./ (b.^2.*cos(ts-t0).^2 + a.^2.*sin(ts-t0).^2));
% xys_els = [rell.*cos(ts) rell.*sin(ts)];
% perim_num = sum(sqrt(sum( (xys_els(1:end-1,:)-xys_els(2:end,:)).^2,2)))

%--------------------------------------------------------------------------
% Construct all the points, and then manually construct the triangulation
% (based on the assumed form)
cylpts = [];
rell   = a*b*sqrt( 1 ./ (b.^2.*cos(ts-t0).^2 + a.^2.*sin(ts-t0).^2));
for n = 1:length(zs)           
    cylpts = [cylpts; rell.*cos(ts) rell.*sin(ts) zs(n)*ones(length(ts),1)];
end
cyltri = manu_tribnd(cylpts);

if plot_flg == 1
    figure
    trisurf(cyltri,cylpts(:,1),cylpts(:,2),cylpts(:,3), ...
        'FaceColor','cyan','FaceAlpha', 1) % ,'linestyle','none');
    axis equal
end


%--------------------------------------------------------------------------
% Define the integrand for the calculation of the length of a curve. For
% polar curves its sqrt( r'^2 + r^2)
% 
% so for an ellipse:
%   r(theta) = a*b*(b^2*cos(theta)^2 + a^2*sin(theta)^2)^-1/2
%
% so r^prime(theta) = -1/2*a*b*(b^2*cos(theta)^2 + a^2*sin(theta)^2)^-3/2 *
%            (-2*b^2*cos(theta)*sin(theta) + 2*a^2*sin(theta)*cos(theta))
%--------------------------------------------------------------------------
function out = polar_integrand(theta,a,b,t0)

rtheta = a*b*sqrt( 1 ./ (b.^2.*cos(theta-t0).^2 + a.^2.*sin(theta-t0).^2));
rprime = -1/2*a*b*(b.^2.*cos(theta-t0).^2 + a.^2.*sin(theta-t0).^2).^(-3/2) .* ...
    (-2*b^2*cos(theta-t0).*sin(theta-t0) + 2*a^2*sin(theta-t0).*cos(theta-t0));
%----------------
out    = sqrt(rtheta.^2 + rprime.^2 );
