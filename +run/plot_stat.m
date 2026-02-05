function plot_stat(NAMES,TENSOR_SZ,TYPE,TIME_ALL,TIME,ERR,DIM,alg,mark,col,TYPE_REF,SAVEIT)

    [NALG,N]    = size(TIME);
    Npoly       = 0;
    Nrat        = 0;
    Nnrat       = 0;
    switch TYPE_REF
        case 'A'
            TYPE_REF_str = 'all';
            idx_keep = 1:N;
        case 'AS'
            TYPE_REF_str = 'all sorted';
            %
            [~,idx_tensor] = sort(TENSOR_SZ,'ascend');
            NAMES       = NAMES(idx_tensor);
            TENSOR_SZ   = TENSOR_SZ(idx_tensor);
            TIME_ALL    = TIME_ALL(:,idx_tensor);
            TIME        = TIME(:,idx_tensor);
            ERR         = ERR(:,idx_tensor);
            DIM         = DIM(:,idx_tensor);
            TYPE        = TYPE(idx_tensor);
            %
            idx_keep = [];
            for ii = 1:N
                if any(strcmp(TYPE{ii},'polynomial'))
                    idx_keep = [idx_keep ii];
                    Npoly    = Npoly+1;
                end
            end
            for ii = 1:N
                if any(strcmp(TYPE{ii},'rational'))
                    idx_keep = [idx_keep ii];
                    Nrat     = Nrat+1;
                end
            end
            for ii = 1:N
                if any(strcmp(TYPE{ii},'irrational'))
                    idx_keep = [idx_keep ii];
                    Nnrat    = Nnrat+1;
                end
            end
        case 'P'
            TYPE_REF_str = 'polynomial';
            idx_keep = [];
            for ii = 1:N
                if any(strcmp(TYPE{ii},'polynomial'))
                    idx_keep = [idx_keep ii];
                end
            end
        case 'R'
            TYPE_REF_str = 'rational';
            idx_keep = [];
            for ii = 1:N
                if any(strcmp(TYPE{ii},'rational'))
                    idx_keep = [idx_keep ii];
                end
            end
        case 'PR'
            TYPE_REF_str = 'polynomial \& rational';
            idx_keep = [];
            for ii = 1:N
                if any(strcmp(TYPE{ii},'polynomial'))
                    idx_keep = [idx_keep ii];
                end
                if any(strcmp(TYPE{ii},'rational'))
                    idx_keep = [idx_keep ii];
                end
            end
        case 'I'
            TYPE_REF_str = 'irrational';
            idx_keep = [];
            for ii = 1:N
                if any(strcmp(TYPE{ii},'irrational'))
                    idx_keep = [idx_keep ii];
                end
            end
    end
    NAMES       = NAMES(idx_keep);
    TENSOR_SZ   = TENSOR_SZ(idx_keep);
    TIME_ALL    = TIME_ALL(:,idx_keep);
    TIME        = TIME(:,idx_keep);
    ERR         = ERR(:,idx_keep);
    DIM         = DIM(:,idx_keep);
    TYPE        = TYPE(:,idx_keep);
    [NALG,N]    = size(TIME);
    %%% TIMES
    % >> size-time ALL
    figure, hold on, grid on, %axis tight
    xval = TENSOR_SZ;
    yval = TIME_ALL;
    ymin = min(yval,[],1);
    for ii = 1:NALG
        plot(xval,yval(ii,:),mark{ii},'MarkerEdgeColor',col(ii,:),'DisplayName',alg{ii})            
        plot_best(ii,xval,yval,ymin,N,mark,col)
        plot_not_cvg(xval,yval(ii,:),min(yval(:)),mark{ii},col(ii,:))
    end
    %legend('show','location','eastoutside')
    legend('show','location','eastoutside','Interpreter','latex');%,'NumColumns',4)
    set(gca,'xtick',[1e-3 1e-2 1 10 100 1e3],'xticklabel',{'1 [KB]' '10 [KB]' '1 [MB]' '10 [MB]' '100 [MB]' '1 [GB]'},'TickLabelInterpreter','latex','XScale','log')
    set(gca,'ytick',[1e-3 1 60 60*15 3600 3600*10],'yticklabel',{'1 [ms]' '1 [s]' '1 [min]' '15 [min]' '1 [h]' '10 [h]'},'TickLabelInterpreter','latex','YScale','log')
    xlabel('\textbf{Tensor size}','interpreter','latex')
    ylabel('\textbf{Time (average over all tests)}','interpreter','latex')
    title(['\textbf{Computational cost vs. tensor size (' TYPE_REF_str ' cases)}'],'interpreter','latex')
    [ranking,ranking_str] = rank_methods(yval);
    drawnow
    if SAVEIT; mdspack.plot.figSavePDF(['tex_pdf/figures/stat/size_time_all_' TYPE_REF],.6); end
    
    % >> size-time BEST
    figure, hold on, grid on, %axis tight
    xval = TENSOR_SZ;
    yval = TIME;
    ymin = min(yval,[],1);
    for ii = 1:NALG
        plot(xval,yval(ii,:),mark{ii},'MarkerEdgeColor',col(ii,:),'DisplayName',alg{ii})
        plot_best(ii,xval,yval,ymin,N,mark,col)
        plot_not_cvg(xval,yval(ii,:),min(yval(:)),mark{ii},col(ii,:))
    end
    %legend('show','location','southoutside','Interpreter','latex','NumColumns',4)
    legend('show','location','eastoutside','Interpreter','latex');%,'NumColumns',4)
    set(gca,'xtick',[1e-3 1e-2 1 10 100 1e3],'xticklabel',{'1 [KB]' '10 [KB]' '1 [MB]' '10 [MB]' '100 [MB]' '1 [GB]'},'TickLabelInterpreter','latex','XScale','log')
    set(gca,'ytick',[1e-3 1 60 60*15 3600 3600*10],'yticklabel',{'1 [ms]' '1 [s]' '1 [min]' '15 [min]' '1 [h]' '10 [h]'},'TickLabelInterpreter','latex','YScale','log')
    xlabel('\textbf{Tensor size}','interpreter','latex')
    ylabel('\textbf{Time (best configuration)}','interpreter','latex')
    title(['\textbf{Computational cost vs. tensor size (' TYPE_REF_str ' cases)}'],'interpreter','latex')
    drawnow
    if SAVEIT; mdspack.plot.figSavePDF(['tex_pdf/figures/stat/size_time_' TYPE_REF],.6); end
    
    % >> Time vs Time
    figure, hold on, grid on, axis tight
    yval = TIME;
    xval = min(yval,[],1);
    ymin = xval;%min(yval,[],1);
    for ii = 1:NALG
        %alg_ratio = sprintf([alg{ii} ' \n ' ranking_str{ii}]);
        %alg_ratio = regexprep(alg_ratio,'%','\%');
        plot(xval,yval(ii,:),mark{ii},'MarkerEdgeColor',col(ii,:),'DisplayName',alg{ii})
        plot_best(ii,xval,yval,ymin,N,mark,col)
        plot_not_cvg(xval,yval(ii,:),min(yval(:)),mark{ii},col(ii,:))
    end
    %plot_area(xval,[1.5 2 3 11])
    %legend('show','location','southoutside','Interpreter','latex','NumColumns',4)
    legend('show','location','eastoutside','Interpreter','latex');%,'NumColumns',4)
    set(gca,'xtick',[1e-3 1 60 60*15 3600 3600*10],'xticklabel',{'1 [ms]' '1 [s]' '1 [min]' '15 [min]' '1 [h]' '10 [h]'},'TickLabelInterpreter','latex','XScale','log')
    set(gca,'ytick',[1e-3 1 60 60*15 3600 3600*10],'yticklabel',{'1 [ms]' '1 [s]' '1 [min]' '15 [min]' '1 [h]' '10 [h]'},'TickLabelInterpreter','latex','YScale','log')
    xlabel('\textbf{Best time}','interpreter','latex')
    ylabel('\textbf{Time (best configuration)}','interpreter','latex')
    title(['\textbf{Time vs. Best time (' TYPE_REF_str ' cases)}'],'interpreter','latex')
    drawnow
    if SAVEIT; mdspack.plot.figSavePDF(['tex_pdf/figures/stat/time_time_' TYPE_REF],.6); end
    
    %%% ERRORS
    % >> size-err BEST
    figure, hold on, grid on, %axis tight
    xval = TENSOR_SZ;
    yval = ERR;
    ymin = min(yval,[],1);
    for ii = 1:NALG
        plot(xval,yval(ii,:),mark{ii},'MarkerEdgeColor',col(ii,:),'DisplayName',alg{ii})
        plot_best(ii,xval,yval,ymin,N,mark,col)
        plot_not_cvg(xval,yval(ii,:),min(yval(:)),mark{ii},col(ii,:))
    end
    %legend('show','location','eastoutside','Interpreter','latex')
    legend('show','location','southoutside','Interpreter','latex','NumColumns',4)
    set(gca,'xtick',[1e-3 1e-2 1 10 100 1e3],'xticklabel',{'1 [KB]' '10 [KB]' '1 [MB]' '10 [MB]' '100 [MB]' '1 [GB]'},'TickLabelInterpreter','latex','XScale','log')
    set(gca,'TickLabelInterpreter','latex','YScale','log')
    xlabel('\textbf{Tensor size}','interpreter','latex')
    ylabel('\textbf{RMS error (best configuration)}','interpreter','latex')
    title(['\textbf{Error vs. tensor size (' TYPE_REF_str ' cases)}'],'interpreter','latex')
    drawnow
    if SAVEIT; mdspack.plot.figSavePDF(['tex_pdf/figures/stat/size_err_' TYPE_REF],.6); end
    
    % >> cases-err BEST
    figure, hold on, grid on, axis tight
    xval = 1:N;
    yval = ERR;
    ymin = min(yval,[],1);
    for ii = 1:NALG
        plot(xval,ERR(ii,:),mark{ii},'MarkerEdgeColor',col(ii,:),'DisplayName',alg{ii})
        plot_best(ii,xval,yval,ymin,N,mark,col)
        plot_not_cvg(xval,yval(ii,:),min(yval(:)),mark{ii},col(ii,:))
    end
    legend('show','location','southoutside','Interpreter','latex','NumColumns',4)
    set(gca,'TickLabelInterpreter','latex','YScale','log')
    xlabel('\textbf{Example}','interpreter','latex')
    ylabel('\textbf{RMS error (best configuration)}','interpreter','latex')
    title(['\textbf{Error vs. \# (' TYPE_REF_str ' cases)}'],'interpreter','latex')
    set(gca,'xtick',1:N,'xticklabel',NAMES,'TickLabelInterpreter','latex')
    drawnow
    if SAVEIT; mdspack.plot.figSavePDF(['tex_pdf/figures/stat/case_err_' TYPE_REF],.6); end
    
    % >> Err vs Err 
    %ERR_mini = min(ERR,[],1);
    figure, hold on, grid on, axis tight
    yval = ERR;
    xval = min(yval,[],1);
    ymin = xval;%min(yval,[],1);
    for ii = 1:NALG
        % ratio       = 100*sum(ERR(ii,:)-ERR_mini==0)/N;
        % alg_ratio   = string([alg{ii} ' - ' num2str(ratio) '\%']);
        plot(xval,yval(ii,:),mark{ii},'MarkerEdgeColor',col(ii,:),'DisplayName',alg{ii})
        plot_best(ii,xval,yval,ymin,N,mark,col)
        plot_not_cvg(xval,yval(ii,:),min(yval(:)),mark{ii},col(ii,:))
    end
    %plot_area(xval,[1.5 2 3 11])
    %legend('show','location','eastoutside','Interpreter','latex')
    legend('show','location','southoutside','Interpreter','latex','NumColumns',4)
    set(gca,'TickLabelInterpreter','latex','XScale','log','YScale','log')
    xlabel('\textbf{Best RMS error (best configuration)}','interpreter','latex')
    ylabel('\textbf{RMS error (best configuration)}','interpreter','latex')
    title(['\textbf{Error vs. Best error (' TYPE_REF_str ' cases)}'],'interpreter','latex')
    drawnow
    if SAVEIT; mdspack.plot.figSavePDF(['tex_pdf/figures/stat/err_err_' TYPE_REF],.6); end

    % >> RADAR
    figure, hold on, axis equal
    fun = @(x) max(1-log10(x./min(x,[],1)),0.1);
    %fun = @(x) max((max(x,[],1)-x)./(max(x,[],1)-min(x,[],1)),0.1);
    plot_radar(fun,TIME_ALL,NAMES,mark,col,alg)
    if strcmp(TYPE_REF,'AS'); plot_bar([Npoly,Nrat,Nnrat],log10(TENSOR_SZ),TYPE); end
    title({['\textbf{Computational cost average of all configurations (' TYPE_REF_str ' cases)}']; func2str(fun)},'interpreter','latex')
    if SAVEIT; axis off; mdspack.plot.figSavePDF(['tex_pdf/figures/stat/radar_time_all_' TYPE_REF],1); end
    %
    figure, hold on, axis equal
    fun = @(x) max(1-log10(x./min(x,[],1)),0.1);
    %fun = @(x) max((max(x,[],1)-x)./(max(x,[],1)-min(x,[],1)),0.1);
    plot_radar(fun,TIME,NAMES,mark,col,alg)
    if strcmp(TYPE_REF,'AS'); plot_bar([Npoly,Nrat,Nnrat],log10(TENSOR_SZ),TYPE); end
    title({['\textbf{Computational cost of the best configuration (' TYPE_REF_str ' cases)}']; func2str(fun)},'interpreter','latex')
    if SAVEIT; axis off; mdspack.plot.figSavePDF(['tex_pdf/figures/stat/radar_time_' TYPE_REF],1); end
    %
    figure, hold on, axis equal    
    fun = @(x) max(1-log10(x./min(x,[],1)),0.1);
    fun = @(x) max(1-log10(x./min(x,[],1))/1e1,0.1);
    %fun = @(x) max((max(x,[],1)-x)./(max(x,[],1)-min(x,[],1)),0.1);
    plot_radar(fun,ERR,NAMES,mark,col,alg)
    if strcmp(TYPE_REF,'AS'); plot_bar([Npoly,Nrat,Nnrat],log10(TENSOR_SZ),TYPE); end
    title({['\textbf{Error of the best configuration (' TYPE_REF_str ' cases)}']; func2str(fun)},'interpreter','latex')
    if SAVEIT; axis off; mdspack.plot.figSavePDF(['tex_pdf/figures/stat/radar_err_' TYPE_REF],1); end
    %
    figure, hold on, axis equal
    %fun = @(x) max(1./x./max(1./x,[],1),.1);
    fun = @(x) max(1-log10(x./min(x,[],1)),0.1);
    plot_radar(fun,DIM,NAMES,mark,col,alg)
    if strcmp(TYPE_REF,'AS'); plot_bar([Npoly,Nrat,Nnrat],log10(TENSOR_SZ),TYPE); end
    title({['\textbf{ROM complexity of the best configuration (' TYPE_REF_str ' cases)}']; func2str(fun)},'interpreter','latex')
    if SAVEIT; axis off; mdspack.plot.figSavePDF(['tex_pdf/figures/stat/radar_dim_' TYPE_REF],1); end

end

%%%
function plot_best(ii,xval,yval,ymin,N,mark,col)

    for jj = 1:N
        if yval(ii,jj) == ymin(jj)
            h=plot(xval(jj),yval(ii,jj),mark{ii},'MarkerEdgeColor',col(ii,:),'MarkerFaceColor',[1 .5 0]);
            h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        end
    end
end

%%% 
function plot_area(xval,fact)    
    fact    = [1; fact(:)];
    col     = parula(numel(fact)+3);
    for ll = 2:numel(fact)
        xx = [min(xval) max(xval)];
        patch([xx(1) xx(end) xx(end) xx(1)], ... 
              [xx(1)*fact(ll-1) xx(end)*fact(ll-1) xx(end)*fact(ll) xx(1)*fact(ll)],col(ll-1),'EdgeColor','none','FaceAlpha',.1,'DisplayName',['$+' num2str(100*(fact(ll)-1)) '\%$']);
        %h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
end

%%% 
function plot_not_cvg(xval,yval,YMIN,mark,col)

    if any(isnan(yval))
        tmp = yval;
        idx = find(isnan(tmp));
        h=plot(xval(idx),YMIN/10*ones(1,numel(idx)),mark,'MarkerEdgeColor',[1 1 1]*.8);%col);
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
end

%%% 
function [ranking,ranking_str] = rank_methods(yval)
    
    [NALG,N]    = size(yval);
    tmp         = yval-min(yval,[],1); 
    [~,idx]     = sort(tmp,1);
    ranking_str = cell(NALG,1);
    for ii = 1:NALG
        for jj = 1:NALG
            ranking(ii,jj)  = 100*sum(idx(jj,:)==ii)/N;
            ranking_str{ii} = [ranking_str{ii} num2str(jj) ': ' num2str(ranking(ii,jj),3) '%% '];
            %ranking_str{ii} = sprintf([ranking_str{ii} num2str(jj) ': ' num2str(ranking(ii,jj),3) ' %% ']);
        end
    end

end

%%% 
function plot_radar(fun,yval,NAMES,mark,col,alg)

    [NALG,N]= size(yval);
    %
    ang_    = linspace(0,2*pi,N+1);
    ang_    = ang_(1:end-1);
    xx_     = real(exp(1i*ang_));
    yy_     = imag(exp(1i*ang_));
    % radius per case
    for ii = 1:N
        h=plot([0.1 1]*xx_(ii),[0.1 1]*yy_(ii),'-','Color',[1 1 1]*.9,'LineWidth',1);
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
    %
    ang_    = linspace(0,2*pi,1e3);
    xx_     = real(exp(1i*ang_));
    yy_     = imag(exp(1i*ang_));
    % 
    %circles = 1+log10([1 .9 .8 .7 .5 .1]);
    scale   = 1;%[1 .9 .8 .7 .5];
    circles = 1+log10(scale);
    %circles = fun(circles)
    for ii = 1:numel(circles)
        h=plot(circles(ii)*xx_,circles(ii)*yy_,'-','Color',[1 1 1]*.9,'LineWidth',1);
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
    %
    ang     = linspace(0,2*pi,N+1);
    ang     = ang(1:end-1);
    xx      = real(exp(1i*ang));
    yy      = imag(exp(1i*ang));
    yval    = fun(yval);
    N       = size(yval,2);
    for ii = 1:NALG
        plot(yval(ii,:).*xx,yval(ii,:).*yy,'-','Color',col(ii,:),'LineWidth',4,'Marker',mark{ii},'MarkerSize',20-2*ii,'DisplayName',alg{ii})
        %plot(yval(ii,:).*xx,yval(ii,:).*yy,'-','Color',col(ii,:),'LineWidth',4,'DisplayName',alg{ii})
    end
    for ii = 1:N
        txt = NAMES{ii};%['\textbf{\#' num2str(idx(ii)) '}'];
        if (ang(ii)<pi/2) || (ang(ii)>3*pi/2)
            text(1.1*xx(ii),1.1*yy(ii),txt,'FontSize',12) 
        else
            text(1.1*xx(ii),1.1*yy(ii),txt,'HorizontalAlignment','right','FontSize',12)
        end
    end
    % % small circles
    % for ii = 1:numel(circles)
    %     if circles(ii) > .1
    %         text(0,circles(ii)-.05,['$' num2str(100*scale(ii)) '\%$'],'Interpreter','latex') 
    %     else
    %         text(0,circles(ii)-.05,['$<' num2str(100*scale(ii)) '\%$'],'Interpreter','latex') 
    %     end
    % end
    legend(alg,'location','southoutside','Interpreter','latex','NumColumns',4)
    set(gca,'TickLabelInterpreter','latex','XLim',[-1 1]*1.3,'YLim',[-1 1]*1.3,'XTickLabel',[],'YTickLabel',[],'xtick',[],'ytick',[])
end

%%% 
function plot_bar(N,TENSOR_SZ,TYPE)
    
    % n       = length(TENSOR_SZ);
    % R       = 1.3;
    % theta   = [linspace(0,2*pi,n)];
    % C       = R*exp(1i*theta);
    % x       = real(C);
    % y       = imag(C);
    % z       = zeros(size(x));
    % col     = TENSOR_SZ;
    % h=surface([x;x],[y;y],[z;z],[col;col],'FaceColor','none','EdgeColor','interp','LineWidth',8);
    % h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    % colorbar
    n       = length(TENSOR_SZ);
    R       = 1.3;
    theta0  = 0;
    idx     = 0;
    for ii = 1:length(N)
        idx     = idx(end) + (1:N(ii));
        theta   = theta0+linspace(0,2*pi*N(ii)/n,N(ii)+1);
        theta0  = theta(end);
        theta   = theta(1:end-1);
        C       = R*exp(1i*theta);
        x       = real(C);
        y       = imag(C);
        z       = zeros(size(x));
        col     = TENSOR_SZ(idx);
        %plot([0 x(end)],[0 y(end)],'k-')
        h=surface([x;x],[y;y],[z;z],[col;col],'FaceColor','none','EdgeColor','interp','LineWidth',8);
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        a=colorbar; 
        a.Label.Interpreter = "latex";
        a.Label.String = '\textbf{Tensor size [Mo] ($10^{x}$-scale)}';
        a.Label.FontSize = 16;
        %
        %text(1.6*real(exp(1i*theta0)),1.6*imag(exp(1i*theta0)),TYPE{idx(end)},'Interpreter','latex') 
    end
end

