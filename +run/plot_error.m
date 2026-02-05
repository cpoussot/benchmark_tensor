function plot_error(NUM,MODULUS,NALG,FIG,infoCas,all,col,Ylimits)

%figure(FIG+ceil(NUM/MODULUS)), hold on, axis tight
ax          = gca();
ax.YGrid    = 'on';
kk          = mod((NUM-1)*NALG,MODULUS*NALG);
% area
ylim(Ylimits); hh = gca;
if mod(NUM,2) == 1
    h=patch([0 NALG NALG 0]+kk+1/2, [min(hh.YLim) min(hh.YLim), max(hh.YLim) max(hh.YLim)],[1 1 1]*.7,'EdgeColor','none','FaceAlpha',.1);
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
% data
for ii = 1:size(all,1)
    kk = kk + 1;
    run.boxplot(kk,all(ii,:),.5,col(ii,:));
end
% Text zone
text(kk-ceil(NALG/2)-1,max(Ylimits)/100,{infoCas.dim_names;infoCas.tensor_names},'FontSize',16,'FontWeight','bold')
%text(kk-ceil(NALG/2)-1,max(Ylimits)/10,infoCas.dim_names,'FontSize',16,'FontWeight','bold')
%text(kk-ceil(NALG/2)-1,max(Ylimits)/100,infoCas.tensor_names,'FontSize',16,'FontWeight','bold')
%
set(gca,'YScale','log')
xlabel('\textbf{Example}','interpreter','latex')
ylabel('\textbf{Abs. error}','interpreter','latex')
title('\textbf{Mismatch abs. error (mean, variance, dispersion)}','interpreter','latex')
