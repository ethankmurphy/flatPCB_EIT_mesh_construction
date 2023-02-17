%--------------------------------------------------------------------------
% 
% get_bound_fine_ell
% 
% This function aims to get a best fitting ellipse to the electrode array,
% and then make it extend a factor fac past any point. 
%
%--------------------------------------------------------------------------
function [a,b,t0] = get_bound_fine_ell(fps,fac,dbg_flg)


%--------------------------------------------------------------------------
% Get the best fitting electrode of the data
[nvec,cent,avec,bvec,svec] = get_nrmal_vec([fps 0*fps(:,1)]);
a0 = svec(1)/svec(2);                % Semi-major axis
b0 = 1;                % Semi-minor axis
t0 = atan2(avec(2),avec(1));    % Rotation to the minor axis
if norm(cent) > 1e-6 
    error('stop: Assume the electrode array is centered')
end

%--------------------------------------------------------------------------
% Get the distance in the a-direction
ts   = atan2(fps(:,2),fps(:,1));
rs   = sqrt(sum(fps.^2,2));
trng = pi/6;
is   =find((abs(ts-t0)<trng) | (abs(ts     -(t0+pi)) < trng) | ...
    (abs(ts+2*pi-t0) < trng) | (abs(ts+2*pi-(t0+pi)) < trng) | ...
    (abs(ts-2*pi-t0) < trng) | (abs(ts-2*pi-(t0+pi)) < trng) );
a    = fac*max(rs(is));
b    = a/a0;

%--------------------------------------------------------------------------
% Ensure the factor is accurate for all the points
rell = a*b*sqrt( 1 ./ (b.^2.*cos(ts-t0).^2 + a.^2.*sin(ts-t0).^2));
facs = rell./rs;
minfac = min(facs);
if minfac < fac %then we need to increase the ellipse a bit
    % We want the minimum factor to be the factor, fac = s*minfac;
    s = fac/minfac;
    a = s*a;
    b = s*b;
    rell = a*b*sqrt( 1 ./ (b.^2.*cos(ts-t0).^2 + a.^2.*sin(ts-t0).^2));
    facs = rell./rs;
    minfac = min(facs);
end

%--------------------------------------------------------------------------
% Debug
if dbg_flg == 1
    disp(['Fac=',num2str(fac),': minimum factor = ',num2str(minfac)])
    svec(1)/svec(2)

    [avec bvec]
    [a0 b0]
    [a b]
    figure;hold on
    plot(fps(:,1),fps(:,2),'.k')
    plot(rell.*cos(ts),rell.*sin(ts),'.r')
    axis equal
end


