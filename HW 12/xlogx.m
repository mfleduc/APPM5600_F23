function res = xlogx(x,varargin)
%% Evaluates xlog(x) for arbitrary base
%   x: Points to evlauate at
%   varargin: One argument: if empty, use natural log. Else, use log base
%   varargin{1}
if isempty(varargin)
    res = x.*log(x);
else
    res = x.*log(x)/log(varargin{1});
end
res(x==0) = 0;

end