function plot_scaled(NUM,infoCas,H,mdl)

figure(NUM),
% Big data
N       = length(mdl);
n       = infoCas.n;
bnd     = infoCas.bound;
%
xlim1   = linspace(min(bnd{1}),max(bnd{1}),52)*.99;
xlim2   = linspace(min(bnd{2}),max(bnd{2}),54)*.99;
imax    = 1;
% if n-2 > 0
%     imax = 1;
% end
[X,Y]   = meshgrid(xlim1,xlim2);
err     = cell(N,1);
tab_red = cell(n,1);
tab_ref = zeros(numel(xlim2),numel(xlim1),imax);
% random point
for i = 1:imax
    if n > 2
        xrnd(i,:) = mlf.rand(n-2,infoCas.bound(3:end),infoCas.domain);
    end
end
% reference
for i = 1:imax
    for ii = 1:numel(xlim1)
        for iii = 1:numel(xlim2)
            x_rnd   = [xlim1(ii) xlim2(iii)];
            if n > 2
                x_rnd   = [x_rnd xrnd(i,:)];
            end
            tab_ref(iii,ii,i)   = H(x_rnd);
        end
    end
end
% approximations
for kk = 1:N
    for i = 1:imax
        for ii = 1:numel(xlim1)
            for iii = 1:numel(xlim2)
                x_rnd   = [xlim1(ii) xlim2(iii)];
                if n > 2
                    x_rnd   = [x_rnd xrnd(i,:)];
                end
                tab_red{kk}(iii,ii,i)   = run.eval(mdl{kk},x_rnd);
            end
        end
    end
    maxTabRef       = 1;%max(abs(tab_ref));
    vec             = abs(tab_ref-tab_red{kk})/maxTabRef;
    vec(isnan(vec)) = 1e12;
    vec(isinf(vec)) = 1e12;
    err{kk}         = log10(vec);
    if any(imag(tab_red{kk})~=0)
        warning('Only real part shown')
        tab_ref     = real(tab_ref);
        tab_red{kk} = real(tab_red{kk});
    end
end

err_min = 1e15;
err_max = -err_min;
for kk = 1:numel(err)
    vec     = err{kk}(:);
    vec(isinf(vec)) = [];
    err_min = min(err_min,min(vec));
    err_max = max(err_max,max(vec));
end

%col = get(gca,'colororder');
for i = 1:size(tab_ref,3)
    kkk     = 0;
    tmp     = tab_ref(:,:,i); tmp = tmp(:);
    zmin    = min(tmp);
    zmax    = max(tmp);
    %dif     = max(abs(zmin),abs(zmax))*5/100;
    %zmin    = zmin-dif;
    %zmax    = zmax+dif;
    clf
    for kk = 1:N
        %
        kkk = kkk + 1;
        subplot(N,2,kkk); hold on, grid on
        surf(X,Y,tab_red{kk}(:,:,i),'EdgeColor','none'), hold on
        surf(X,Y,tab_ref(:,:,i),'EdgeColor','k','FaceColor','none')
        xlabel('$x_1$','Interpreter','latex')
        ylabel('$x_2$','Interpreter','latex')
        title(mdl{kk}.name,'Interpreter','latex')
        set(gca,'TickLabelInterpreter','latex')
        axis tight, zlim([zmin zmax]), view(-30,40)
        %
        kkk = kkk + 1;
        subplot(N,2,kkk); hold on, grid on, axis tight
        imagesc(err{kk}(:,:,i),'XData',xlim1,'YData',xlim2)
        xlim([min(xlim1) max(xlim1)])
        ylim([min(xlim2) max(xlim2)])
        xlabel('$x_1$','Interpreter','latex')
        ylabel('$x_2$','Interpreter','latex')
        title(['{\bf log}(abs. err.)'],'Interpreter','latex')
        set(gca,'TickLabelInterpreter','latex')
        colorbar, clim([err_min err_max])
        %sgtitle([infoCas.name],'FontSize',20,'interpreter','latex')
        if n > 2
            xrnd_str = [];
            for ii = 1:numel(xrnd)
                xrnd_str = [xrnd_str num2str(xrnd(ii)) ' ; '];
            end
            xrnd_str = xrnd_str(1:end-3);
        end
        if n == 3
            sgtitle(['$x_{3}=[' xrnd_str ']$'],'FontSize',20,'interpreter','latex')
        elseif n > 3
            sgtitle(['$x_{3...' num2str(n) '}=[' xrnd_str ']$'],'FontSize',20,'interpreter','latex')
        end
        drawnow
    end
end
