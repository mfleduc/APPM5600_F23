function [root, errMsg, x] = newtonsMethod(f, fp, x0, tol, maxIters, varargin)
%% Code for running Newton's method for rootfinding on a 1D function
%Inputs:
%   f: Function handle, the function we are finding the root of
%   fp: Function handle, the derivative of f
%   x0: Scalar, initial gruess
%   tol: scalar, relative error tolerance
%   maxIters: Scalar, maximum number of iterations allowed
%   varargin: Should be either the multiplicity of the root, or empty.
%Outputs:
%   root: scalar, approximate root location
%   errMsg: Flag that contains information on success/failure of the method
%       0: Success
%       1: Failure, took too many iterations
%       2: Failure, ran off to infinity (i.e. because the function was too
%       flat near the initial guess)
%   x: vector x_n, the sequence of approximations of the root   
if ~isempty(varargin)
     p = varargin{1};
else
    p = 1;
end
x = zeros(1, maxIters);x(1)=x0;
cnt = 1;
del = tol+1;
while del>tol && cnt<maxIters && ~isnan(del) && f(x(cnt))~=0
    dx = p*f(x(cnt))/fp(x(cnt));
    x(cnt+1) = x(cnt)-dx;
    del = abs(dx/(x(cnt+1)));
    cnt=cnt+1;
end
x = x(1:cnt);
root = x(end);
if cnt == maxIters
    errMsg = 1;
elseif isnan(del)
    errMsg = 2;
else
    errMsg = 0;
end
end