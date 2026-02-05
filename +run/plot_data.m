function plot_data(NUM,MODULUS,NALG,FIG,infoCas,data,col,Ylimits)

%figure(FIG+ceil(NUM/MODULUS)), hold on, axis tight
ax          = gca();
ax.YGrid    = 'on';
kk          = mod((NUM-1)*NALG,MODULUS*NALG);
% area
ylim(Ylimits); hh = gca;
if mod(NUM,2) == 1
    h=patch([0 NALG NALG 0]+kk+1/2,10.^[-5 -5, 10 10],[1 1 1]*.7,'EdgeColor','none','FaceAlpha',.1);
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
% data
for jj = 1:size(data,1)
    kk = kk + 1;
    if ~isnan(data(jj))
        bar(kk,data(jj),'facecolor','flat','CData',col(jj,:),'FaceAlpha',.7,'EdgeColor',col(jj,:),'BarWidth',.5);
    else
        plot(kk, max(hh.YLim(:)), '*', 'color',col(jj,:));
        text(kk, min(hh.YLim(:)*10),'NOT CONVERGED','Rotation',90,'Color',col(jj,:),'FontWeight','bold');
    end
end
% Text zone
text(kk-ceil(NALG/2)-1,max(Ylimits)/5,{infoCas.dim_names;infoCas.tensor_names},'FontSize',16,'FontWeight','bold')
%text(kk-ceil(NALG/2)-1,max(Ylimits)/5,infoCas.dim_names,'FontSize',16,'FontWeight','bold')
%text(kk-ceil(NALG/2)-1,max(Ylimits)/10,infoCas.tensor_names,'FontSize',16,'FontWeight','bold')
%
set(gca,'YScale','log')
xlabel('\textbf{Example}','interpreter','latex')
ylabel('\textbf{Complexity (\# parameters)}','interpreter','latex')
title('\textbf{Model complexity}','interpreter','latex')
