function functionList = make_latex_main(CAS,infoCas)

if infoCas.tab_MB < 1;     tensor_size = num2str(infoCas.tab_MB*2^10,3); unit = '\textbf{KB}'; end
if infoCas.tab_MB >= 1;    tensor_size = num2str(infoCas.tab_MB,3);      unit = '\textbf{MB}'; end
if infoCas.tab_MB >= 2^10; tensor_size = num2str(infoCas.tab_MB/2^10,3); unit = '\textbf{GB}'; end

functionList    = [];
functionList    = [functionList '\newpage '];
if strcmp(infoCas.name(1),'$')
    functionList    = [functionList ['\subsection{Function \#' num2str(CAS) ' (${\ord=' num2str(infoCas.n) '}$ variables, tensor size: ' tensor_size ' ' unit ')} $$' infoCas.name(2:end-1) '$$ ' ]];
else
    functionList    = [functionList ['\subsection{Function \#' num2str(CAS) ' (${\ord=' num2str(infoCas.n) '}$ variables, tensor size: ' tensor_size ' ' unit ')} $$\texttt{' infoCas.name '}$$ ' ]];
end
functionList    = [functionList '\subsubsection{Setup and results overview}'];
functionList    = [functionList '\begin{itemize}'];
functionList    = [functionList ['\item Reference: ' infoCas.ref ', ' infoCas.cite]];
if strcmp(infoCas.domain,'R')
    functionList    = [functionList ['\item Domain: $\mathbb{R}$' ]];
elseif strcmp(infoCas.domain,'C')
    functionList    = [functionList ['\item Domain: $\mathbb{C}$' ]];
end

functionList    = [functionList ['\item Tensor size: ' tensor_size ' ' unit ' ($' num2str(infoCas.nip*2) '^{' num2str(infoCas.n) '}$ points)']];
functionList    = [functionList ['\item Bounds: $ ' ]];
for ii = 1:numel(infoCas.bound)
    functionList    = [functionList latex(sym(infoCas.bound{ii})) ];
    if ii < numel(infoCas.bound)
        functionList    = [functionList ' \times '];
    end
end
functionList    = [functionList '$ \end{itemize} '];
functionList    = [functionList ['\begin{table}[H] \centering \input{figures/case_' num2str(CAS) '/table_main.tex} \caption{Function \#' num2str(CAS) ': best model configuration and performances per methods.} \end{table}' ]];
functionList    = [functionList ['\begin{figure}[H] \centering  \includegraphics[width=\textwidth]{figures/case_' num2str(CAS) '/all_stat.pdf} \caption{Function \#' num2str(CAS) ': graphical view of the best model performances.} \end{figure}' ]];
functionList    = [functionList ['\begin{figure}[H] \centering  \includegraphics[width=\textwidth]{figures/case_' num2str(CAS) '/eval_scaled.pdf} \caption{Function \#' num2str(CAS) ': left side, evaluation of the original (mesh) vs. approximated (coloured surface) and right side, absolute errors (in log-scale).} \end{figure}' ]];
% image_latex     = ['figures/stat_3D_CAS' num2str(CAS) '_scaled.pdf' ];
% if exist(['tex_pdf/' image_latex])
%    functionList = [functionList ['\begin{figure}[H] \begin{center}\includegraphics[width=\textwidth]{' image_latex '} \end{center}\caption{Function \#' num2str(CAS) ': left side, evaluation of the original (mesh) vs. approximated (coloured surface) and right side, absolute errors (in log-scale).} \end{figure}' ]];
% end
% 
functionList    = [functionList ['\subsubsection{mLF detailed informations (M1)} \input{figures/case_' num2str(CAS) '/text_loe.tex}']];

