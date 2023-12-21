function [x,errEst, cnt, MiNnorm, errMsg] = stationaryIters( A,b, x0, tol, maxIters, method, varargin )
%Stationary iterations using Jacobi, Gauss-Seidel, or successive
%over-relaxation methods to solve Ax=b for x
%Inputs:
%   A: The matrix, Ax=b
%   b: The y-value. Ax=b
%   x0: Initial guess. 
%   tol: Relative error tolerance
%   maxIters: Max number of iterations allowed
%   method: Which algorithm to use. One of 'jacobi', 'gauss-seidel', or
%       'sor'
%   varargin: If using SOR, the value of \omega
%Outputs:
%   x: Estimate of the solution
%   errEst: Error estimate, c/(1-c)*norm(dx)
%   cnt: Number of iterations
%   MiNnorm: operator norm of M^{-1}N
%   errMsg: 0 if successful, 1 if ran out of iterations, 2 if returned NaN
[n,m] = size(A,[1,2]);
if n~=m
   error( 'A is not square' ) 
end
if strcmpi(method, 'jacobi')%% Pick an algorithm
   M = diag(A).*eye(n);
   L = tril(A)-M;
   U = triu(A)-M;
   N = -(L+U);
elseif strcmpi(method, 'gauss-seidel')
   M = tril(A);
   N = M-A;
elseif strcmpi(method, 'SOR')
    if ~isempty(varargin{1})
        w = varargin{1};
    else
        w = 1.6735;
    end
    D = diag(A).*eye(n);
    L = tril(A)-D;
    U = triu(A)-D;
    M = D+w*L;
    N = (1-w)*D-w*U;
    b = w*b;
end
Mi = inv(M);
MiNnorm = norm(Mi*N, 2);
del = tol+1;
cnt = 0;
while cnt<maxIters && del>tol %Implementing the iteration
    x1 = Mi*(N*x0+b);
    del = norm(x1-x0)/norm(x0);
    cnt=cnt+1;
    xn1 = x0;
    x0 = x1;
end
errMsg = 0;
x=x0;
eigVals = eig( Mi*N );%% Getting teh error estimate
spectralRad = max(abs(eigVals));
errEst = abs(spectralRad/(1-spectralRad))*norm( x0-xn1 );
if cnt<maxIters && all(isfinite(x))
    errMsg = 0;
elseif all(isfinite(x))
    errMsg = 1;
else
    errMsg = 2;
end
end