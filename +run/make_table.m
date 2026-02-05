function Tk = make_table(NAMES,METHOD_LIST,par,dims,time,err_rms)

    NALG        = length(dims);
    name_NUM    = cell(NALG,1); name_NUM{1} = NAMES;
    name_ALG    = METHOD_LIST;
    name_PAR    = fun_texBold(par);%par;%regexprep(num2str(par),'\s*',',');
    name_DIMS   = fun_texBold(dims);
    name_TIME   = fun_texBold(time);
    name_RMSE   = fun_texBold(err_rms);
    %name_MINE   = fun_texBold(err_min);
    %name_MAXE   = fun_texBold(err_max);
    Tk          = table(name_NUM,name_ALG(:),name_PAR,name_DIMS,name_TIME,name_RMSE);%,name_MINE,name_MAXE);
    Tk.Properties.VariableNames = {'$\#$','Alg.','Parameters','Dim.','CPU [s]','RMSE'};%,'min err.','max err.'};

end

%%%
function value_out = fun_texBold(value)

    if isa(value,'cell')
        for ii = 1:numel(value)
            value_out{ii,1} = regexprep(num2str(value{ii}),'\s*',',');
            value_out{ii,1} = regexprep(value_out{ii,1}, 'e\+?(-?\d+)', ' \\cdot 10^{$1}');
            value_out{ii,1} = ['$' value_out{ii,1} '$'];
        end
    else
        [~,idx] = min(value);
        for ii = 1:numel(value)
            if ii == idx
                value_out{ii,1} = ['\mathbf{' num2str(value(ii),2) '}'];
            else
                value_out{ii,1} = [num2str(value(ii),2)];
            end
            value_out{ii,1} = regexprep(value_out{ii,1}, 'e\+?(-?\d+)', ' \\cdot 10^{$1}');
            value_out{ii,1} = ['$' value_out{ii,1} '$'];
        end
    end
    
    % [~,idx]         = min(value);
    % value_out       = num2str(value(ii),2);
    % value_out(idx)  = ['$\mathbf{' num2str(value(idx),2) '}$'];
end