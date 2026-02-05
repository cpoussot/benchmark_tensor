function [alg,par,dims,time,idx,err_rms,err_min,err_max,err_all,time_all] = errors(varargin)

% Big data
N       = nargin-2;
H       = varargin{1};
p_rnd   = varargin{2};
nn      = size(p_rnd,1);
% reference
for ii = 1:nn
    tab_ref(ii,1) = H(p_rnd(ii,:));
end
% approximations
alg     = cell(N,1);
par     = cell(N,1);
err_rms = zeros(N,1);
err_min = zeros(N,1);
err_max = zeros(N,1);
dims    = zeros(N,1);
time    = zeros(N,1);
for kk = 1:N
    alg{kk,1}   = varargin{kk+2}.name;
    par{kk,1}   = varargin{kk+2}.par;
    dims(kk,1)  = varargin{kk+2}.data_sz;
    time(kk,1)  = varargin{kk+2}.time;
    tab_app = zeros(nn,1);
    for ii = 1:nn
        tic;
        tab_app(ii)     = run.eval(varargin{kk+2},p_rnd(ii,:));
        time_all(kk,ii) = toc;
    end
    %
    %err_all(kk,:)   = (abs(tab_ref-tab_app)).';
    err_all(kk,:)   = (abs(tab_ref-tab_app)).'./max(abs(tab_ref(:)));
    err_rms(kk,1)   = rmse(tab_ref,tab_app);
    err_max(kk,1)   = max(err_all(kk,:));
    err_min(kk,1)   = min(err_all(kk,:));
end
[~,idx] = min(err_rms);