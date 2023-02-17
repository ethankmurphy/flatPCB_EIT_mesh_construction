%-------------------------------------------------------------------------------
%
%
%
%-------------------------------------------------------------------------------
function plot_msh_elecs_only(msh,nde_scl,plt_elnum,alpha)

if nargin == 1
    nde_scl   = 1;
    plt_elnum = 0;
elseif nargin == 2
    plt_elnum = 0;
end

% figure
hold on
nds = msh.node*nde_scl;
ecents = get_el_cents(msh)*nde_scl;
for k = 1:length(msh.elec)
    ftmp   = msh.face(msh.elec{k},:);
    fcs    = 1/3*(nds(ftmp(:,1),:)+nds(ftmp(:,2),:)+nds(ftmp(:,3),:));
    ecent  = mean(fcs,1);
    ndstmp = nds;
    trisurf(ftmp,nds(:,1),nds(:,2),nds(:,3), ...
        ... % 'FaceColor','red','FaceAlpha', alpha);
        'FaceColor','red','FaceAlpha', alpha,'edgecolor','red','linestyle','none');
    if plt_elnum == 1
        text(ecents(k,1),ecents(k,2),ecents(k,3)+0.001,num2str(k),'fontsize',12)
    end
end
view(2)
