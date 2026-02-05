function tab = eval(model,p)

tab = NaN;
switch lower(model.method)
    case 'loewner'
        try 
            tab = model.G(p);
        catch e
            fprintf(1,'The identifier was: \n %s\n',e.identifier)
        end
    case 'mds'
        try
            p_str = regexprep(num2str(p,36),'\s*',',');
            eval(['tab = model.G(' p_str ');'])
        catch e
            fprintf(1,'The identifier was: \n %s\n',e.identifier)
        end
    case 'paaa'
        try
            tab = model.G(p);
        catch e
            fprintf(1,'The identifier was: \n %s\n',e.identifier)
        end
    case 'kan1'
        try 
            method  = model.data{1};
            xmin    = model.data{2};
            xmax    = model.data{3};
            ymin    = model.data{4};
            ymax    = model.data{5};
            fnB     = model.data{6};
            fnT     = model.data{7};
            if ~(any(isnan(fnB(:))) | any(isnan(fnT(:))))
                if method < 4
                    tab = modelKA_basisC(p, xmin, xmax, ymin, ymax, fnB, fnT );
                else
                    tab = modelKA_linear(p, xmin, xmax, ymin, ymax, fnB, fnT );
                end
            end
        catch e
            fprintf(1,'The identifier was: \n %s\n',e.identifier)
        end
    case 'tf'
        W1      = model.data{1};
        b1      = model.data{2};
        W2      = model.data{3};
        b2      = model.data{4};
        relu    = @(x) (x+abs(x))./2; %max(0,x)
        p       = p(:).';
        tab     = relu(p*W1+b1)*W2+b2;
    %otherwise 
    %    error('Unknown model structure')
end
