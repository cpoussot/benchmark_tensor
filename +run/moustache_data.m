function [miny, maxy, meany, medy, iQ1, iQ3] = moustache_data(y)
miny    = min(y);
maxy    = max(y);
meany   = mean(y);
medy    = median(y);

N       = length(y);
iQ1     = ceil(N/4);
iQ3     = ceil(3*N/4);
end