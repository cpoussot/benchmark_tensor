clear variables; close all; clc; format short
%%% Addpath of all the methods
add_path_third_parties

%%% Chose the directory for result save
RESULT_PATH     = '/Users/charles/Library/CloudStorage/ProtonDrive-charles.poussot@proton.me-folder/Research/Benchmarks/mLF_evaluation/vtest/results';

%%% Chose method list & set tensor size limit
METHOD_LIST     = {'mlf1' 'mlf2' 'mdspack' 'kan1' 'paaa' 'paaalr'}; 
TENSOR_MAX      = [ inf    inf    inf       30     30     30]; % in MB

%%% Parameter combinations for each methods
% MLF1 (https://github.com/cpoussot/mLF)
paramMLF1.tol   = [.5 1e-1 1e-2 1e-3 1e-4 1e-6 1e-9 1e-10 1e-11 1e-12 1e-13 1e-14];
paramMLF1.null  = {'svd0' 'qr0' 'mldivide0'};
paramMLF1.comb  = mlf.combinations_dim([length(paramMLF1.tol) length(paramMLF1.null)]);
% MLF2 (https://github.com/cpoussot/mLF)
paramMLF2.tol   = 1e-15;
paramMLF2.null  = paramMLF1.null;
paramMLF2.comb  = mlf.combinations_dim([length(paramMLF2.tol) length(paramMLF2.null)]);
% MDS (https://mordigitalsystems.fr/)
paramMDS.method = {'R'};
paramMDS.tolk   = [1e-2 1e-4 1e-6 1e-8 1e-10 1e-12 1e-14 1e-15 -1];
paramMDS.tol    = [0 paramMLF1.tol];
paramMDS.comb   = mlf.combinations_dim([length(paramMDS.method) length(paramMDS.tolk) length(paramMDS.tol)]);
% KAN (https://github.com/andrewpolar)
paramKAN.method = 1:4;
paramKAN.alpha  = [.95 1];  % damping factor for iterative parameter update (also called learning rate) 
paramKAN.lambda = .01;      % Tikhonov regularisation parameter for Gauss-Newton method 
paramKAN.Nrun   = 50;       % num. of runs through data
paramKAN.n      = [4 6 10]; % num. of nodes bottom
paramKAN.q      = [4 6 12]; % num. of nodes top
paramKAN.comb   = mlf.combinations_dim([length(paramKAN.method) length(paramKAN.alpha) length(paramKAN.Nrun) length(paramKAN.lambda) length(paramKAN.n) length(paramKAN.q)]);
% pAAA (https://github.com/lbalicki/parametric-AAA)
paramAAA.tol    = [1e-3 1e-6 1e-9];
paramAAA.comb   = mlf.combinations_dim([length(paramAAA.tol)]);
% pAAA-LR (https://github.com/lbalicki/parametric-AAA)
paramAAALR.tol  = paramAAA.tol;
paramAAALR.rank = [2 3 4 5];
paramAAALR.comb = mlf.combinations_dim([length(paramAAALR.tol) length(paramAAALR.tol)]);
% TensorFlow (https://www.tensorflow.org/?hl=fr)
% Nota vailable in this version yet

%%% Number of random draw and examples number
NTEST           = 500;
spaceCAS        = 2:3;%1:50;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% NO NEED TO CHANGE FROM HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for CAS = spaceCAS%1:50 
    %%% Load model and extract variables
    [H,infoCas] = mlf.examples(CAS)
    n           = infoCas.n;
    p_c         = infoCas.p_c;
    p_r         = infoCas.p_r;
    %%% Build tensor grid
    tic
    [y,x,dim]   = mlf.make_tab_vec(H,p_c,p_r);
    tab         = mlf.vec2mat(y,dim);
    N           = length(y);
    toc, pause(1)
    %%% MLF (A1) [A/G/P-V, 2025]
    if any(strcmpi(METHOD_LIST,'mlf1')) && (infoCas.tab_MB <= TENSOR_MAX(1))
        if ~exist([RESULT_PATH '/mlf1']); mkdir([RESULT_PATH '/mlf1']); end
        clear mdl; arg_ = [];
        for ii = 1:size(paramMLF1.comb,1) 
            tic;
            MLF1_tol            = paramMLF1.tol(paramMLF1.comb(ii,1));
            MLF1_null           = paramMLF1.null(paramMLF1.comb(ii,2));
            %
            opt                 = [];
            opt.ord_tol         = MLF1_tol;
            opt.method_null     = MLF1_null{1};
            opt.method          = 'rec';
            opt.ord_obj         = [];
            opt.ord_show        = false;
            opt.data_min        = true;
            [G,info]            = mlf.alg1(tab,p_c,p_r,opt);
            %
            mdl{ii}.name        = 'A/G/P-V 2025 (A1)';
            mdl{ii}.method      = 'loewner';
            mdl{ii}.par         = [MLF1_tol paramMLF1.comb(ii,2)];
            mdl{ii}.G           = G;
            mdl{ii}.data        = {info.pc info.w info.c};
            mdl{ii}.data_sz     = numel(ones(prod(info.ord+1),2+n));
            mdl{ii}.info        = info;
            mdl{ii}.time        = toc; 
            arg_                = [arg_ 'mdl{' num2str(ii) '},'] ;
            %
            arg = arg_(1:end-1);
            save([RESULT_PATH '/mlf1/cas_' num2str(CAS) '_mlf1'],'infoCas','H','mdl','arg')
        end
    end
    %%% MLF (A2) [A/G/P-V, 2025]
    if any(strcmpi(METHOD_LIST,'mlf2')) && (infoCas.tab_MB <= TENSOR_MAX(2))
        if ~exist([RESULT_PATH '/mlf2']); mkdir([RESULT_PATH '/mlf2']); end
        clear mdl; arg_ = [];
        MLF2_maxIter    = 20;
        for ii = 1:size(paramMLF2.comb,1)
            tic;
            MLF2_tol            = paramMLF2.tol(paramMLF2.comb(ii,1));
            MLF2_null           = paramMLF2.null(paramMLF2.comb(ii,2));
            %
            opt                 = [];
            opt.tol             = MLF2_tol;
            opt.method_null     = MLF2_null{1};
            opt.max_iter        = MLF2_maxIter;
            opt.method          = 'rec';
            try 
                [G,info]        = mlf.alg2(tab,p_c,p_r,opt);
                data            = {info.pc info.w info.c};
                data_sz         = numel(ones(prod(info.ord+1),2+n));
                time            = toc;
            catch e
                fprintf(1,'The identifier was: \n %s\n',e.identifier)
                toc;
                G               = NaN;
                data            = NaN;
                data_sz         = NaN;
                time            = NaN;
                info            = NaN;
            end
            mdl{ii}.name        = 'A/G/P-V 2025 (A2)';
            mdl{ii}.method      = 'loewner';
            mdl{ii}.par         = [MLF2_tol paramMLF2.comb(ii,2)];
            mdl{ii}.G           = G;
            mdl{ii}.data        = data;
            mdl{ii}.data_sz     = data_sz;
            mdl{ii}.time        = time;
            arg_                = [arg_ 'mdl{' num2str(ii) '},'] ;
            %
            arg = arg_(1:end-1);
            save([RESULT_PATH '/mlf2/cas_' num2str(CAS) '_mlf2'],'infoCas','H','mdl','arg')
        end
    end
    %%% MDSPACK [https://mordigitalsystems.fr/]
    if any(strcmpi(METHOD_LIST,'mdspack')) && (infoCas.tab_MB <= TENSOR_MAX(3))
        if ~exist([RESULT_PATH '/mdspack']); mkdir([RESULT_PATH '/mdspack']); end
        clear mdl; arg_ = [];
        for ii = 1:size(paramMDS.comb,1)
            clear ip
            tic;
            MDS_method          = paramMDS.method(paramMDS.comb(ii,1));
            MDS_tolk            = paramMDS.tolk(paramMDS.comb(ii,2));
            MDS_tol             = paramMDS.tol(paramMDS.comb(ii,3));
            % 
            opt                 = [];
            opt.method          = MDS_method{1};
            opt.tolk            = MDS_tolk;
            if MDS_tol == 0
                for iii = 1:n; ip{iii} = [p_c{iii} p_r{iii}]; end
                tabr = tab;
            else
                ord             = mlf.compute_order(p_c,p_r,tab,MDS_tol,[],5,false);
                [pc,pr,~,~,tabr]= mlf.points_selection(p_c,p_r,tab,ord,true);
                for iii = 1:n; ip{iii,1} = [pc{iii}(:); pr{iii}(:)]; end
            end
            [G,info]            = mdspack.ratapp(tabr,ip,opt);
            %
            mdl{ii}.name        = 'MDSPACK v1.1.0';
            mdl{ii}.method      = 'mds';
            mdl{ii}.par         = [MDS_tolk paramMDS.comb(ii,2)];
            mdl{ii}.G           = G;
            mdl{ii}.data        = info.r;
            mdl{ii}.data_sz     = numel(info.r);
            mdl{ii}.time        = toc;
            arg_                = [arg_ 'mdl{' num2str(ii) '},'] ;
            %
            arg = arg_(1:end-1);
            save([RESULT_PATH '/mdspack/cas_' num2str(CAS) '_mdspack'],'infoCas','H','mdl','arg')
        end
    end
    %%% KAN1 [P/P, 2025]
    if any(strcmpi(METHOD_LIST,'kan1')) && (infoCas.tab_MB <= TENSOR_MAX(4))
        if ~exist([RESULT_PATH '/kan1']); mkdir([RESULT_PATH '/kan1']); end
        clear mdl; arg_ = [];
        xmin    = min(infoCas.xmin);
        xmax    = max(infoCas.xmax);
        ymin    = min(y(:)); 
        ymax    = max(y(:));
        for ii = 1:size(paramKAN.comb,1)
            tic;
            KAN_method  = paramKAN.method(paramKAN.comb(ii,1));
            KAN_alpha   = paramKAN.alpha(paramKAN.comb(ii,2));
            KAN_Nrun    = paramKAN.Nrun(paramKAN.comb(ii,3));
            KAN_lambda  = paramKAN.lambda(paramKAN.comb(ii,4));
            KAN_n       = paramKAN.n(paramKAN.comb(ii,5));
            KAN_q       = paramKAN.q(paramKAN.comb(ii,6));
            KAN_p       = 2*KAN_n+1;
            Nid         = floor(N/2)+1;
            identID     = 1; lab = ones(N,1);
            verifID     = 2; lab(Nid:end) = 2;
            try 
                [fnB0,fnT0] = buildKA_init( n, KAN_n, KAN_q, KAN_p, ymin, ymax );
                if (KAN_method == 1) %. basis functions - cubic splines, identification method - Gauss-Newton
                    [yhat_all, fnB, fnT, RMSE, t_min_all, t_max_all ] = solveMinGauss( x, y, lab, identID, verifID, KAN_alpha, KAN_lambda, KAN_Nrun, xmin, xmax, ymin, ymax, fnB0, fnT0 );
                elseif (KAN_method == 2) %. basis functions - cubic splines, identification method - Newton-Kaczmarz, standard
                    [yhat_all, fnB, fnT, RMSE, t_min_all, t_max_all ] = buildKA_basisC(x, y, lab, identID, verifID, KAN_alpha, KAN_Nrun, xmin, xmax, ymin, ymax, fnB0, fnT0 );
                elseif (KAN_method == 3) %. basis functions - cubic splines, identification method - Newton-Kaczmarz, accelerated
                    [yhat_all, fnB, fnT, RMSE, t_min_all, t_max_all ] = buildKA_basisA(x, y, lab, identID, verifID, KAN_alpha, KAN_Nrun, xmin, xmax, ymin, ymax, fnB0, fnT0 );
                elseif (KAN_method == 4) %. basis functions - piecewise-linear, identification method - Newton-Kaczmarz, standard
                    [yhat_all, fnB, fnT, RMSE, t_min_all, t_max_all ] = buildKA_linear(x, y, lab, identID, verifID, KAN_alpha, KAN_Nrun, xmin, xmax, ymin, ymax, fnB0, fnT0 );
                end
                data_sz = 4+numel(fnB)+numel(fnT);
                time    = toc;
            catch e
                fprintf(1,'The identifier was: \n %s\n',e.identifier)
                toc;
                data_sz = NaN;
                fnB     = NaN;
                fnT     = NaN;
                time    = NaN;
            end
            mdl{ii}.name    = 'P/P 2025';
            mdl{ii}.method  = 'kan1';
            mdl{ii}.par     = [KAN_method KAN_alpha KAN_Nrun KAN_lambda KAN_n KAN_q KAN_p];
            mdl{ii}.G       = [];
            mdl{ii}.data    = {KAN_method xmin xmax ymin ymax fnB fnT};
            mdl{ii}.data_sz = data_sz;
            mdl{ii}.time    = time;
            arg_            = [arg_ 'mdl{' num2str(ii) '},'] ;
            %
            arg = arg_(1:end-1);
            save([RESULT_PATH '/kan1/cas_' num2str(CAS) '_kan1'],'infoCas','H','mdl','arg')
        end
    end
    %%% pAAA [B/G, 2025]
    if any(strcmpi(METHOD_LIST,'paaa')) && (infoCas.tab_MB <= TENSOR_MAX(5))
        if ~exist([RESULT_PATH '/paaa']); mkdir([RESULT_PATH '/paaa']); end
        clear mdl; arg_ = [];
        for ii = 1:size(paramAAA.comb,1)
            clear ip
            tic;
            AAA_maxIter     = 20;
            opt.max_iter    = AAA_maxIter;
            %
            AAA_tol         = paramAAA.tol(paramAAA.comb(ii,1));
            for iii = 1:infoCas.n; ip{iii} = [p_c{iii} p_r{iii}]; end
            try 
                [G,info]    = paaa(tab,ip,AAA_tol,opt);
                Gfun        = @(p) G.eval(p);
                data_sz     = prod(G.orders+1)*(n+2);
                time        = toc;
            catch e
                fprintf(1,'The identifier was: \n %s\n',e.identifier)
                toc;
                Gfun        = NaN;
                G           = NaN;
                data_sz     = NaN;
                time        = NaN;
                info        = NaN;
            end
            mdl{ii}.name    = 'C-R/B/G 2023';
            mdl{ii}.method  = 'paaa';
            mdl{ii}.par     = [AAA_tol AAA_maxIter];
            mdl{ii}.G       = Gfun;
            mdl{ii}.data    = G;
            mdl{ii}.data_sz = data_sz;
            mdl{ii}.time    = time;
            arg_            = [arg_ 'mdl{' num2str(ii) '},'] ;
            %
            arg = arg_(1:end-1);
            save([RESULT_PATH '/paaa/cas_' num2str(CAS) '_paaa'],'infoCas','H','mdl','arg')
        end
    end
    %%% pAAA-LR [B/G, 2025]
    if any(strcmpi(METHOD_LIST,'paaalr')) && (infoCas.tab_MB <= TENSOR_MAX(6))
        if ~exist([RESULT_PATH '/paaalr']); mkdir([RESULT_PATH '/paaalr']); end
        clear mdl; arg_ = [];
        for ii = 1:size(paramAAALR.comb,1)
            clear ip
            tic;
            AAALR_maxIter   = 20;
            opt.max_iter    = AAALR_maxIter;
            %
            AAALR_tol       = paramAAALR.tol(paramAAALR.comb(ii,1));
            AAALR_rank      = paramAAALR.rank(paramAAALR.comb(ii,2));
            for iii = 1:infoCas.n; ip{iii} = [p_c{iii} p_r{iii}]; end
            try 
                [G,info]    = lr_paaa(tab,ip,AAALR_tol,AAALR_rank,opt);
                Gfun        = @(p) G.eval(p);
                data_sz     = prod(G.orders+1)*(n+2);
                time        = toc;
            catch e
                fprintf(1,'The identifier was: \n %s\n',e.identifier)
                toc;
                Gfun        = NaN;
                G           = NaN;
                data_sz     = NaN;
                time        = NaN;
                info        = NaN;
            end
            mdl{ii}.name    = 'B/G 2025';
            mdl{ii}.method  = 'paaa';
            mdl{ii}.par     = [AAALR_tol AAALR_maxIter AAALR_rank];
            mdl{ii}.G       = Gfun;
            mdl{ii}.data    = G;
            mdl{ii}.data_sz = data_sz;
            mdl{ii}.time    = time;
            arg_            = [arg_ 'mdl{' num2str(ii) '},'] ;
            %
            arg = arg_(1:end-1);
            save([RESULT_PATH '/paaalr/cas_' num2str(CAS) '_paaalr'],'infoCas','H','mdl','arg')
        end
    end
    % %%% TensorFlow
    % if any(strcmp(lower(METHOD_LIST),'tensorflow'))
    %     clear mdl; arg_ = [];
    %     %name_py = "run_tensorflow/run2D.py";
    %     %res = pyrunfile(name_py,"z",x=3,y=2)
    %     % for ii = 1%:size(paramTF.comb,1)
    %     % 
    %     % end
    % end
end
