clear variables; close all; clc; format short
%%% Addpath of all the methods
add_path_third_parties

%%% Chose the directory for result save 
% /!\ Same as "RESULT_PATH", line 6 in "start_compare_step1.m"
RESULT_PATH = '/Users/charles/Library/CloudStorage/ProtonDrive-charles.poussot@proton.me-folder/Research/Benchmarks/mLF_evaluation/vtest/results';
TEX_PATH    = '/Users/charles/Library/CloudStorage/ProtonDrive-charles.poussot@proton.me-folder/Research/Benchmarks/mLF_evaluation/vtest/tex_pdf';

%%% Chose method list 
METHOD_LIST = {'mlf1' 'mlf2' 'mdspack' 'kan1' 'paaa' 'paaalr'}; 

%%% Variables
spaceCas    = 1:50;
MODULUS     = 10;
PLOT_2D     = true;
PLOT_ERR    = true;
PLOT_TIME   = true;
PLOT_DATA   = true;
PLOT_ALL    = true;
PLOT_STAT   = true;
MAKE_TEX    = true;
SAVEIT      = true;

%%% Colors & shapes
col         = run.load_colors();
mark        = {'s' 'd' 'o' '^' '>' 'v' '<'};
%
NALG        = numel(METHOD_LIST);
FIG_ERR     = 1e3;
FIG_TIME    = 1e4;
FIG_DATA    = 1e5;
FIG_ALL     = 1e7;
fun_time    = @(x) mean( x(x>0 & ~isnan(x)));
%%% Constant random seed to ensure reproducibility (at least on a given MATLAB setting)
rng(1712)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% NO NEED TO CHANGE FROM HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kk  = 0;
if ~exist(TEX_PATH); mkdir(TEX_PATH); end
for CAS = spaceCas
    %
    if ~exist([TEX_PATH '/figures/case_' num2str(CAS)]); mkdir([TEX_PATH '/figures/case_' num2str(CAS)]), end 
    if ~exist([TEX_PATH '/figures/stat']); mkdir([TEX_PATH '/figures/stat']), end 
    %
    kk = kk + 1;
    %%% Load model 
    [H,infoCas] = mlf.examples(CAS)
    NAMES{kk}   = infoCas.cas_names;
    TYPE{kk}    = infoCas.tag{1};
    if MAKE_TEX
        tex_main    = run.make_latex_main(CAS,infoCas);
        writelines(tex_main,[TEX_PATH '/figures/case_' num2str(CAS) '/text_main.tex'])
        %
        tex_slide = ['\begin{center}' infoCas.name '\end{center} \begin{itemize}' ...  
                     '\item Reference: ' infoCas.ref ...
                     '\item Tensor size: ' infoCas.tensor_names '\end{itemize}'];
        writelines(tex_slide,[TEX_PATH '/figures/case_' num2str(CAS) '/text_slide.tex'])
    end

    %%% Load results
    for ii = 1:NALG
        fileName        = [RESULT_PATH '/' METHOD_LIST{ii} '/cas_' num2str(CAS) '_' METHOD_LIST{ii} '.mat'];
        fileBestName    = [RESULT_PATH '/' METHOD_LIST{ii} '/cas_' num2str(CAS) '_' METHOD_LIST{ii} '_best.mat'];
        if exist(fileBestName)
            tmp             = load(fileBestName);
            tmp             = tmp.best;
            idx_best        = tmp.idx;
            time_all{ii,:}  = tmp.time;
            %
            alg{ii,1}       = tmp.alg{idx_best};
            par{ii,1}       = tmp.par{idx_best};
            %
            dims(ii,1)      = tmp.dims(idx_best);
            time(ii,1)      = tmp.time(idx_best);
            err_rms(ii,1)   = tmp.err_rms(idx_best);
            err_min(ii,1)   = tmp.err_min(idx_best);
            err_max(ii,1)   = tmp.err_max(idx_best);
            err(ii,:)       = tmp.err_all(idx_best,:);
            time_eval(ii,:) = tmp.time_all(idx_best,:);
            %
            tmp             = load(fileName);
            mdl_opt{ii}     = tmp.mdl{idx_best};
        else
            time_all{ii,1}  = NaN;
            par{ii,1}       = NaN;
            %
            dims(ii,1)      = NaN;
            time(ii,1)      = NaN;
            err_rms(ii,1)   = NaN;
            err_min(ii,1)   = NaN;
            err_max(ii,1)   = NaN;
            err(ii,:)       = NaN;
            time_eval(ii,:) = NaN;
            %
            mdl_opt{ii}         = []; 
            mdl_opt{ii}.method  = '';
            mdl_opt{ii}.name    = run.default_name(METHOD_LIST{ii});
        end
    end
    for ii = 1:NALG
        if isempty(alg{ii,1}); alg{ii,1}   = ''; end
    end
    %%% Collect data
    TENSOR_SZ(kk)   = infoCas.tab_MB;
    for ii = 1:NALG
        TIME_ALL(ii,kk)     = fun_time(time_all{ii});
        TIME_BEST(ii,kk)    = time(ii);
        DIM_BEST(ii,kk)     = dims(ii);
        ERR_BEST(ii,kk)     = err_rms(ii);%fun_time(err_rms(ii));
    end

    %%% Mini table
    if MAKE_TEX
        Tk = run.make_table(NAMES,alg,par,dims,time,err_rms);
        mlf.table2latex_tabular(Tk,[TEX_PATH '/figures/case_' num2str(CAS) '/table_main.tex']);
        tex_loe = run.make_latex_loe(mdl_opt{1});
        writelines(tex_loe,[TEX_PATH '/figures/case_' num2str(CAS) '/text_loe.tex'])
    end

    %%% 2D Plot best only
    if PLOT_2D
        run.plot_scaled(CAS,infoCas,H,mdl_opt);
        if SAVEIT; mlf.figSavePDF([TEX_PATH '/figures/case_' num2str(CAS) '/eval_scaled'],1.4); end
    end
    
    %%% Plot errors
    if PLOT_ERR 
        figure(FIG_ERR+ceil(kk/MODULUS)), hold on, axis tight
        run.plot_error(kk,MODULUS,NALG,FIG_ERR,infoCas,err,col,[1e-19 1e5])
        legend(alg,'Interpreter','latex','Location','southoutside','NumColumns',4)
        if (mod(kk,MODULUS) == 0)
            set(gca,'xtick',(1:NALG:NALG*MODULUS)+1,'xticklabel',NAMES(kk-MODULUS+1:kk),'TickLabelInterpreter','latex')
            if SAVEIT; mlf.figSavePDF([TEX_PATH '/figures/stat/err_' num2str(kk/MODULUS)],.5); end
        elseif (kk == length(spaceCas))
            set(gca,'xtick',(1:NALG:NALG*MODULUS)+1,'xticklabel',NAMES(kk-mod(kk,MODULUS)+1:kk),'TickLabelInterpreter','latex')
            if SAVEIT; mlf.figSavePDF([TEX_PATH '/figures/stat/err_' num2str(ceil(kk/MODULUS))],.5); end
        end
        drawnow, pause(.5)
    end

    %%% Plot time
    if PLOT_TIME 
        figure(FIG_TIME+ceil(kk/MODULUS)), hold on, axis tight
        run.plot_time(kk,MODULUS,NALG,FIG_TIME,infoCas,time,col,10.^[-4 6])
        legend(alg,'Interpreter','latex','Location','southoutside','NumColumns',4)
        if (mod(kk,MODULUS) == 0)
            set(gca,'xtick',(1:NALG:NALG*MODULUS)+1,'xticklabel',NAMES(kk-MODULUS+1:kk),'TickLabelInterpreter','latex')
            if SAVEIT; mlf.figSavePDF([TEX_PATH '/figures/stat/time_' num2str(kk/MODULUS)],.5); end
        elseif (kk == length(spaceCas))
            set(gca,'xtick',(1:NALG:NALG*MODULUS)+1,'xticklabel',NAMES(kk-mod(kk,MODULUS)+1:kk),'TickLabelInterpreter','latex')
            if SAVEIT; mlf.figSavePDF([TEX_PATH '/figures/stat/time_' num2str(ceil(kk/MODULUS))],.5); end
        end
        drawnow, pause(.5)
    end

    %%% Plot data
    if PLOT_DATA
        figure(FIG_DATA+ceil(kk/MODULUS)), hold on, axis tight
        run.plot_data(kk,MODULUS,NALG,FIG_DATA,infoCas,dims,col,[1 1e6])
        legend(alg,'Interpreter','latex','Location','southoutside','NumColumns',4)
        if (mod(kk,MODULUS) == 0)
            set(gca,'xtick',(1:NALG:NALG*MODULUS)+1,'xticklabel',NAMES(kk-MODULUS+1:kk),'TickLabelInterpreter','latex')
            if SAVEIT; mlf.figSavePDF([TEX_PATH '/figures/stat/data_' num2str(kk/MODULUS)],.5); end
        elseif (kk == length(spaceCas))
            set(gca,'xtick',(1:NALG:NALG*MODULUS)+1,'xticklabel',NAMES(kk-mod(kk,MODULUS)+1:kk),'TickLabelInterpreter','latex')
            if SAVEIT; mlf.figSavePDF([TEX_PATH '/figures/stat/data_' num2str(ceil(kk/MODULUS))],.5); end
        end
        drawnow, pause(.5)
    end

    %%% Plot all
    if PLOT_ALL
        figure(FIG_ALL+kk)
        subplot(131), grid on, hold on
        run.plot_error(0,1,NALG,FIG_ERR,infoCas,err,col,[1e-19 1e5])
        xlabel('','interpreter','latex'), set(gca,'xtick',[],'TickLabelInterpreter','latex')
        title('\textbf{Mismatch abs. error}')
        subplot(132), grid on, hold on
        run.plot_time(0,1,NALG,FIG_TIME,infoCas,time,col,10.^[-4 6])
        xlabel('','interpreter','latex'), set(gca,'xtick',[],'TickLabelInterpreter','latex')
        subplot(133), grid on, hold on
        run.plot_data(0,1,NALG,FIG_DATA,infoCas,dims,col,[1 1e6])
        xlabel('','interpreter','latex'), set(gca,'xtick',[],'TickLabelInterpreter','latex')
        drawnow, pause(.5)
        legend(alg,'Interpreter','latex','Location','North')
        if SAVEIT; mlf.figSavePDF([TEX_PATH '/figures/case_' num2str(CAS) '/all_stat'],.5); end
    end
end

%%
% %%% Plot size-time
% if PLOT_STAT
%     % >> 
%     TYPE_REF = 'A';
%     %TYPE_REF = 'AS';
%     %TYPE_REF = 'R';
%     %TYPE_REF = 'P';
%     %TYPE_REF = 'I';
%     %TYPE_REF = 'PR';
%     run.plot_stat(NAMES,TENSOR_SZ,TYPE,TIME_ALL,TIME_BEST,ERR_BEST,DIM_BEST,alg,mark,col,TYPE_REF,SAVEIT)
% end