set(groot,'DefaultFigurePosition', [300 100 1000 700]);
set(groot,'defaultlinelinewidth',2.5)
set(groot,'defaultlinemarkersize',14)
set(groot,'defaultaxesfontsize',18)
set(groot,'defaulttextinterpreter','latex')
set(groot,'defaultlegendinterpreter','latex')

% MLF
addpath('/Users/charles/Documents/GIT/mlf')
% MDSPACK
mdspack_ver = 'v1.1.0';
setenv('MDSHOME','/Users/charles/Documents/MDS')
addpath(['/Users/charles/Documents/MDS/mdspack/MDSPACK/osx/' mdspack_ver '/API/matlab/'])
addpath(['/Users/charles/Documents/MDS/mdspack/MDSPACK/osx/' mdspack_ver '/bin'])
% KAN
addpath('/Users/charles/Documents/GIT/_others/kan-polar')
% pAAA & TT
addpath('/Users/charles/Documents/GIT/_others/tensor_toolbox')
addpath('/Users/charles/Documents/GIT/_others/parametric-AAA')
addpath('/Users/charles/Documents/GIT/_others/parametric-AAA/barycentric_forms')
addpath('/Users/charles/Documents/GIT/_others/parametric-AAA/lr_paaa')
addpath('/Users/charles/Documents/GIT/_others/parametric-AAA/paaa')
addpath('/Users/charles/Documents/GIT/_others/parametric-AAA/sv_paaa')
addpath('/Users/charles/Documents/GIT/_others/parametric-AAA/utils')
addpath('/Users/charles/Documents/GIT/_others/parametric-AAA/utils/cauchy_mats')
addpath('/Users/charles/Documents/GIT/_others/parametric-AAA/utils/realness') 
