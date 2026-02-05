clear variables; close all; clc; format short
%%% Addpath of all the methods
add_path_third_parties

%%% Chose the directory for result save 
% /!\ Same as "RESULT_PATH", line 6 in "start_compare_step1.m"
RESULT_PATH     = '/Users/charles/Library/CloudStorage/ProtonDrive-charles.poussot@proton.me-folder/Research/Benchmarks/mLF_evaluation/vtest/results';

%%% Chose method list 
METHOD_LIST     = {'mlf1' 'mlf2' 'mdspack' 'kan1' 'paaa' 'paaalr'}; 

%%% Number of random draw and examples number
NTEST           = 500;
spaceCAS        = 1:50;

%%% Constant random seed to ensure reproducibility (at least on a given MATLAB setting)
rng(1712)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% NO NEED TO CHANGE FROM HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for CAS = spaceCAS
    %%% Load model 
    [H,infoCas] = mlf.examples(CAS)
    %%% Random draw
    clear p_rnd
    for ii = 1:NTEST
        p_rnd(ii,:) = mlf.rand(infoCas.n,infoCas.bound,infoCas.domain);
    end
    %%% Errors by methods and save best configurations
    for ii = 1:numel(METHOD_LIST)
        fileName    = [RESULT_PATH '/' METHOD_LIST{ii} '/cas_' num2str(CAS) '_' METHOD_LIST{ii} '.mat'];
        if exist(fileName)
            fileName
            data            = load(fileName);
            mdl             = data.mdl;
            eval(['[alg,par,dims,time,idx,err_rms,err_min,err_max,err_all,time_all] = run.errors(H,p_rnd,' data.arg ');'])
            % algo
            best.alg        = alg;
            best.par        = par;
            best.dims       = dims;
            best.time       = time;
            % evaluation
            best.idx        = idx;
            best.err_rms    = err_rms;
            best.err_min    = err_min;
            best.err_max    = err_max;
            best.err_all    = err_all;
            best.time_all   = time_all;
            save([RESULT_PATH '/' METHOD_LIST{ii} '/cas_' num2str(CAS) '_' METHOD_LIST{ii} '_best'],'best');
        end
    end
end
