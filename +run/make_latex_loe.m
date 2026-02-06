%function latexList = make_latex_loe(p_c,c,w,Var,Lag,Bary)
function latexList = make_latex_loe(mdl_opt)

FUN_SYM = @(x) vpa(x,4);
p_c     = mdl_opt.data{1};
w       = mdl_opt.data{2};
c       = mdl_opt.data{3};
lag     = mdl_opt.info.lag;
p_r     = mdl_opt.info.pr;
n       = numel(p_c);
N       = numel(c);
Nmax    = 30;%20
%%% Generic stuff
% vars    = '(';
% ord     = [];
% for ii = 1:numel(p_c)
%     vars    = [vars 's_{' num2str(ii) '},'];
%     ord     = [ord num2str(length(p_c{ii})-1) ','];
%     k_i(ii) = length(p_c{ii});
% end
% vars    = [vars(1:end-1) ')'];
% ord     = ord(1:end-1);
vars    = '(';
vars_z  = '(';
vars_y  = '(';
for ii = 1:n
    vars    = [vars   's_{' num2str(ii) '},'];
    vars_z  = [vars_z 'z_{' num2str(ii) '},'];
    vars_y  = [vars_y 'y_{' num2str(ii) '},'];
    k_i(ii) = length(p_c{ii});
    %q_i(ii) = length(p_r{ii});
end
%N_i     = size(tab);
vars    = [vars(1:end-1)   ')'];
vars_z  = [vars_z(1:end-1) ')'];
vars_y  = [vars_y(1:end-1) ')'];
%
latexList   = [];
if N <= Nmax
    %%% KST
    [Var,Lag,Bary]  = mlf.decoupling(p_c,lag);
    d = 1;
    for ii = 1:numel(p_c)%Var
        d = d.*Lag{ii}.*Var{ii}; % den
    end

    %%% Interpolation points
    % right
    latexList = [latexList '\noindent \textbf{Right interpolation points} ($k_l=' latex(sym(k_i)) '$, where $l=1,\cdots,\ord$):' ];
    latexList = [latexList '$$ \begin{array}{rcl}'];
    for ii = 1:length(p_c)
        latexList = [latexList ['\lan{' num2str(ii) '} &=& ' latex(sym(p_c{ii})) '\\' ]];
        %latexList = [latexList ['\mun{' num2str(ii) '} &=& ' latex(sym(p_r{ii})) '\\' ]];
    end
    latexList = [latexList '\end{array} $$' ];
    
    %%% Loewner c, w and c.w
    latexList   = [latexList '\noindent \textbf{Lagrangian weights}: '];
    top         = {'rc' 'rw' 'rcw'};
    MAT         = sym([c w w.*c]);
    for ii = 1:3
        tmp_data(:,ii) = [top{ii}; MAT(:,ii)];
    end
    latexList   = [latexList [ '$$' latex(FUN_SYM(tmp_data)) '$$' ]];

    %%% Transfer function in Lagrangian
    [glag,ilag] = mlf.tf_lagrangian(p_c,w,c,false);
    [num,den]   = numden(simplify(glag));
    %num         = ilag.num_coeff;
    %den         = ilag.den_coeff;
    %num         = simplify(sum(num.*ilag.basis));
    %den         = simplify(sum(den.*ilag.basis));
    d0 = den; for iii = 1:n; d0 = subs(d0,['s' num2str(iii)],0); end; num = num/d0; den = den/d0;
    latexList   = [latexList '\noindent \textbf{Lagrangian form} '];
    latexList   = [latexList '(basis, numerator and denominator coefficients):'];
    latexList   = [latexList '$$\left(\begin{array}{ccc}\mathcal{B}_\textrm{lag}' vars ' & \bN_\textrm{lag} &\bD_\textrm{lag}\end{array}\right) =$$ $$' latex(FUN_SYM([ilag.basis ilag.num_coeff ilag.den_coeff])) '.$$' ];
    latexList   = [latexList '\noindent The corresponding function is:'];
    latexList   = [latexList '$$\begin{array}{rcl}\bG_{\textrm{lag}}' vars ' &=& \dfrac{\bn_{\textrm{lag}}' vars '}{\bd_{\textrm{lag}}' vars '}\\ && \\&=& \dfrac{\sum_{\textrm{row}} \bN_\textrm{lag} \odot\mathcal{B}^{-1}_\textrm{lag}' vars '}{\sum_{\textrm{row}} \bD_\textrm{lag} \odot\mathcal{B}^{-1}_\textrm{lag}' vars '}, \end{array}$$'];
    latexList   = [latexList '\noindent where,\\'];
    latexList   = [latexList '$\bn_{\textrm{lag}}' vars ' = ' latex(FUN_SYM(num)) '$ \\~~\\'];
    latexList   = [latexList '$\bd_{\textrm{lag}}' vars ' = ' latex(FUN_SYM(den)) '$ \\~~\\'];

    %%% Transfer function in Monomial
    [gmon,imon] = mlf.tf_monomial(p_c,w,c,false);
    %[num,den]   = numden(gmon);
    num         = imon.num_coeff;
    den         = imon.den_coeff;
    num         = simplify(sum(num.*imon.basis));
    den         = simplify(sum(den.*imon.basis));
    d0 = den; for iii = 1:n; d0 = subs(d0,['s' num2str(iii)],0); end; num = num/d0; den = den/d0;
    latexList   = [latexList '\noindent \textbf{Monomial form}'];
    latexList   = [latexList ' (basis, numerator and denominator coefficients - entries $<10^{-12}$ removed):'];
    latexList   = [latexList '$$\left(\begin{array}{ccc}\mathcal{B}_\textrm{mon}' vars ' & \bN_\textrm{mon} &\bD_\textrm{mon}\end{array}\right) =$$ $$' latex(FUN_SYM([imon.basis imon.num_coeff imon.den_coeff])) '$$' ];
    latexList   = [latexList '\noindent The corresponding function is:'];
    latexList   = [latexList '$$\begin{array}{rcl}\bG_{\textrm{mon}}' vars ' &=& \dfrac{\bn_{\textrm{mon}}' vars '}{\bd_{\textrm{mon}}' vars '}\\ && \\&=& \dfrac{\sum_{\textrm{row}} \bN_\textrm{mon} \odot \mathcal{B}_\textrm{mon}' vars '}{\sum_{\textrm{row}} \bD_\textrm{mon} \odot\mathcal{B}_\textrm{mon}' vars '},  \end{array}$$'];
    latexList   = [latexList '\noindent where,\\'];
    latexList   = [latexList '$\bn_{\textrm{mon}}' vars ' = ' latex(FUN_SYM(num)) '$ \\~~\\'];
    latexList   = [latexList '$\bd_{\textrm{mon}}' vars ' = ' latex(FUN_SYM(den)) '$ \\~~\\'];

    %%% KST decoupling
    [Bary,Lag,Cx]  = mlf.decoupling(p_c,lag);
    eval(['maxLength = length(Cx.d' num2str(n) ');'])
    if maxLength > 10
        latexList = [latexList '\begin{landscape} '];
    end
    latexList = [latexList '\noindent \textbf{KST equivalent decoupling pattern}'];
    % si's
    latexList = [latexList ' (Barycentric weights $\bc^{\var{l}}$): $$\begin{array}{rclll}'];
    str_k     = [];
    tmp_k     = [];
    for ii = numel(p_c):-1:1
        eval(['tmp = Cx.s' num2str(ii) ';']);
        if ii == numel(p_c)
            latexList   = [latexList ['\var{' num2str(ii) '}&: & ' latex(FUN_SYM(tmp)) '& \textrm{vec}(.) & := \textbf{Bary}(\var{' num2str(ii) '}) \\' ]];
        else
            latexList   = [latexList ['\var{' num2str(ii) '}&: & ' latex(FUN_SYM(tmp)) '& \textrm{vec}(.) \otimes ' str_k ' & := \textbf{Bary}(\var{' num2str(ii) '}) \\' ]];
        end
        tmp_k   = [tmp_k ['k_{' num2str(ii) '}']];
        str_k   = ['\bone_{' tmp_k '}'];
    end
    latexList = [latexList ['\end{array}$$' ]];
    % % w's
    % latexList = [latexList ['~\\ Data weights $\bw^{\var{l}}$: $$\begin{array}{rcll}']];
    % str_k     = [];
    % tmp_k     = [];
    % for ii = numel(p_c):-1:1
    %     eval(['tmp = Cx.w' num2str(ii) ';']);
    %     if ii == numel(p_c)
    %         latexList   = [latexList ['\var{' num2str(ii) '}&: & ' latex(FUN_SYM(tmp)) '& \textrm{vec}(.) \\' ]];
    %     else
    %         latexList   = [latexList ['\var{' num2str(ii) '}&: & ' latex(FUN_SYM(tmp)) '& \textrm{vec}(.) \otimes ' str_k '\\' ]];
    %     end
    %     tmp_k   = [tmp_k ['k_{' num2str(ii) '}']];
    %     str_k   = ['\bone_{' tmp_k '}'];
    % end
    % latexList = [latexList '\end{array}$$' ];
    % % scalings
    % latexList = [latexList '~\\ Normalization scalings $d^{\var{l}}$: $$\begin{array}{rcll}'];
    % str_k     = [];
    % tmp_k     = [];
    % for ii = numel(p_c):-1:1
    %     eval(['tmp = Cx.d' num2str(ii) ';']);
    %     if ii == numel(p_c)
    %         latexList   = [latexList ['\var{' num2str(ii) '}&: & ' latex(FUN_SYM(tmp)) '& \textrm{vec}(.) \\' ]];
    %     else
    %         latexList   = [latexList ['\var{' num2str(ii) '}&: & ' latex(FUN_SYM(tmp)) '& \textrm{vec}(.) \otimes ' str_k '\\' ]];
    %     end
    %     tmp_k   = [tmp_k ['k_{' num2str(ii) '}']];
    %     str_k   = ['\bone_{' tmp_k '}'];
    % end
    % latexList = [latexList '\end{array}$$' ];
    if maxLength > 10
        latexList = [latexList '\end{landscape} '];
    end
    % 
    latexList = [latexList '~\\ Then, with the above notations, one defines the following univariate vector functions:  $$ \left\{ \begin{array}{rcl}'];
    for ii = 1:n
        latexList = [latexList ['\bPhi_{' num2str(ii) '}(\var{' num2str(ii) '}) &:=& \textbf{Bary}(\var{' num2str(ii) '}) \odot \mathbf{Lag}(\var{' num2str(ii) '}) \\']];
    end
    latexList = [latexList '\end{array} \right. $$'];

    %%% Resulting KST
    latexList   = [latexList '\noindent The corresponding function is:'];
    latexList   = [latexList '$$\begin{array}{rcl}\bG_{\textrm{kst}}' vars ' &=& \dfrac{\bn_{\textrm{kst}}' vars '}{\bd_{\textrm{kst}}' vars '}\\ && \\ &=& \dfrac{\sum_{\text{rows}} \bw \odot \bPhi_{1}(\var{1}) \odot \cdots \odot\bPhi_{' num2str(n) '}(\var{' num2str(n) '})}{\sum_{\text{rows}} \bPhi_{1}(\var{1}) \odot \cdots \odot\bPhi_{' num2str(n) '}(\var{' num2str(n) '})} . \end{array}$$'];

    %%% KST univariate functions 
    [gkst1,kst] = mlf.make_singleVarFun(p_c,p_r,Cx);
    gkst1       = kst.f;
    latexList   = [latexList '~\\ \noindent \textbf{KST-like univariate functions} '];
    latexList   = [latexList '(equivalent scaled univariate functions $\bphi_{1,\cdots,' num2str(length(gkst1)) '}$): '];
    latexListF  = '$$\left\{\begin{array}{rcrcl}';
    idx         = [];
    %inc         = 1;
    for ii = 1:length(gkst1)
        eq_ii   = latex(FUN_SYM(simplify(gkst1{ii})));
        ZVAR    = ['z_{' num2str(ii) '} &=&'];
        if length(eq_ii) < 100
            %latexListF  = [latexListF [ZVAR '\bphi_{' num2str(ii) '}(s_{' num2str(vars_1(ii)) '}) &=& ' eq_ii '\\' ]];
            latexListF  = [latexListF [ZVAR '\bphi_{' num2str(ii) '}(s_{' num2str(ii) '}) &=& ' eq_ii '\\' ]];
        else
            idx         = [idx ii];
            [num,den]   = numden(simplify(gkst1{ii}));
            d0 = den; for iii = 1:n; d0 = subs(d0,['s' num2str(iii)],0); end; num = num/d0; den = den/d0;
            nn{ii}      = latex(FUN_SYM(simplify(num)));
            dd{ii}      = latex(FUN_SYM(simplify(den)));
            %latexListF  = [latexListF [ZVAR '\bphi_{' num2str(ii) '}(s_{' num2str(vars_1(ii)) '}) &=& \frac{\bn_{' num2str(ii) '}}{\bd_{' num2str(ii) '}} \\' ]];
            latexListF  = [latexListF [ZVAR '\bphi_{' num2str(ii) '}(s_{' num2str(ii) '}) &=& \frac{\bn_{' num2str(ii) '}}{\bd_{' num2str(ii) '}} \\' ]];
        end
    end
    latexListF  = [latexListF '\end{array} \right. .$$' ];
    latexListF  = strrep(latexListF,'\frac','\cfrac');
    latexList   = [latexList latexListF];
    if ~isempty(idx)
        latexList   = [latexList '\noindent where, \\ '];
        for ii = 1:length(idx)
            latexList   = [latexList '\noindent $\bn_{' num2str(idx(ii)) '}=' nn{idx(ii)} '$ and \\ '];
            latexList   = [latexList '\noindent $\bd_{' num2str(idx(ii)) '}=' dd{idx(ii)} '$, \\ '];
        end
    end

else
    %%% Interpolation points
    latexList = [latexList '\noindent \textbf{Right interpolation points}: $k_l=' latex(sym(k_i)) '$, where $l=1,\cdots,\ord$.' ];
    latexList = [latexList '$$ \begin{array}{rcl}'];
    for ii = 1:length(p_c)
        latexList = [latexList ['\lan{' num2str(ii) '} &\in& \IC^{' num2str(numel(p_c{ii})) '} \text{ , linearly spaced between bounds}\\' ]];
        %latexList = [latexList ['\mun{' num2str(ii) '} &=& ' latex(sym(p_r{ii})) '\\' ]];
    end
    latexList   = [latexList '\end{array} $$' ];
    latexList   = [latexList '\noindent \textbf{$\ord$-D Loewner matrix, barycentric weights and Lagrangian basis}:'];
    latexList   = [latexList '$$ \begin{array}{rcl}'];
    n           = num2str(length(c));
    latexList   = [latexList ['\IL & \in & \IC^{' n ' \times ' n '}\\' ]];
    latexList   = [latexList ['\bc & \in & \IC^{' n '}\\' ]];
    latexList   = [latexList ['\bw & \in & \IC^{' n '}\\' ]];
    latexList   = [latexList ['\bc\odot \bw & \in & \IC^{' n '}\\' ]];
    latexList   = [latexList ['\mathbf{Lag}' vars ' & \in & \IC^{' n '}\\' ]];
    latexList   = [latexList '\end{array} $$' ];
end

%%% Some minor text-style
latexList = strrep(latexList,'\mathrm{rc}','\bc');
latexList = strrep(latexList,'\mathrm{rw}','\bw');
latexList = strrep(latexList,'\mathrm{rcw}','\bc\odot\bw');
latexList = strrep(latexList,'\mathrm{lag}','\mathbf{Lag}^{-1}');
% % Replace bary with \bc 
% for ii = 1:numel(p_c)
%     latexList = strrep(latexList,['\mathrm{bary}_{' num2str(ii) '}'],['\bc^{\var{' num2str(ii) '}}\odot \mathbf{Lag}(\var{' num2str(ii) '})']);
% end
% Replace s with \var 
for ii = 1:numel(p_c)
    latexList = strrep(latexList,['s_{' num2str(ii) '}'],['\var{' num2str(ii) '}']);
end
% 10^
latexList = regexprep(latexList, 'e\+?(-?\d+)', ' \\cdot 10^{$1}');
% \frac
% latexList = strrep(latexList,'\frac','\cfrac');

%%% Equivalent NN
if N < 20
    %
    latexList   = [latexList '\noindent \textbf{Connection with Neural Networks} (equivalent numerator $\bn_{\textrm{lag}}$ representation):\\ '];
    latexNN     = mlf.make_latex_NN_lag(p_c,double(w.*c)); 
    latexList   = [latexList ['\begin{figure}[H]\begin{center} \scalebox{.7}{' latexNN '} \caption{Equivalent NN representation of the numerator $\bn_{\textrm{lag}}$.}\end{center}\end{figure}']];
    % %
    % latexList   = [latexList '\noindent \textbf{Connection with Neural Networks} (equivalent denominator $\bd_{\textrm{lag}}$ representation):\\ '];
    % latexNN     = mlf.make_latex_NN_lag(p_c,double(c)); 
    % latexList   = [latexList ['\begin{figure}[H]\begin{center} \scalebox{.7}{' latexNN '} \caption{Equivalent NN representation of the denominator $\bd_{\textrm{lag}}$.}\end{center}\end{figure}']];
end
