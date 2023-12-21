function [x0,ierr, resid, iterates] = ConjGrad(A,b, tol, maxIters)
%% Conjugate gradient
%Inputs
%   A: An SPD matrix, solving Ax=b
%   b: RHS of Ax=b
%   tol: Relative error tolerance
%   maxIters: Maximum number of iterations allowed
%Outputs
%   x0: Approximate solution
%   ierr: Error message
%       0: Successful
%       1: Failed, too many iterations
%       Note: This will only return 1 if you select fewer iterations than
%       size(A,1) and the method does not converge in the desired
%       number of iterations
%       2: A NaN showed up
%   resid: Relative residuals at each step
%   iterates: VAlues of x_k at each step

iterates = zeros(length(b), maxIters);
x0 = zeros(size(b(:)));
p = b;
% p = W*tmp;
r = p;
cnt = 0;
dx = 1+tol;
maxIters = min(maxIters, max(size(A)));
resid = zeros(1, maxIters);

while dx > tol && cnt<maxIters
    cnt=cnt+1;
    Ap = A*p;
    alpha = norm(r)^2/sum(p'*Ap);
    x1 = x0+alpha*p;
    r1 = r-alpha*Ap;
    beta = (norm(r1)/norm(r))^2;
    p1 = beta*p+r1;
    dx = norm(x0-x1)/norm(x1);
    resid(cnt) = dx;
    %% Update the stuff
    x0 = x1;
    r = r1;
    p = p1;
    iterates(:,cnt) = x1;
end

if cnt == maxIters && maxIters<max(size(A))
    ierr = 1 ;
    fprintf('Conjugate gradient: Desired tolerance not reached, loop reached max iterations\n');
elseif isnan(dx)
    ierr = 2;
    fprintf('Conjugate gradient: Desired tolerance not reached, NaN encountered\n');
else
    ierr = 0;
    fprintf('Conjugate gradient: Reached desired tolerance in %d iterations\n',cnt);
end
