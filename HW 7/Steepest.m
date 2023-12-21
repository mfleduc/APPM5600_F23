function [x,ierr,resid,iterates] = Steepest(A,b,x0,tol,maxIters)
%% Steepest descent method 
%Inputs
%   A: An SPD matrix, solving Ax=b
%   b: RHS of Ax=b
%   x0: Initial guess
%   tol: Relative error tolerance
%   maxIters: Maximum number of iterations allowed
%Outputs
%   x: Approximate solution
%   ierr: Error message
%       0: Successful
%       1: Failed, too many iterations
%       2: NaN made an appearance somewhere it shouldn't have
%   resid: Relative residuals at each step
%   iterates: Values of x_k at each step
iterates = zeros(length(b), maxIters);
del = tol+1;
cnt = 0;
resid = zeros(1, maxIters);
while cnt<maxIters && del>tol && ~isnan(del)
    cnt = cnt+1;
    r = b-A*x0;
    alpha = norm(r)^2/(r'*A*r);
    x = x0+alpha*r;
    del = norm(x-x0)/norm(x);
    resid(cnt) = del;
    x0 = x;
    iterates(:,cnt) = x;
end
if cnt==maxIters 
    ierr = 1;
    fprintf('Steepest Descent: Desired tolerance not reached, loop reached max iterations\n');
elseif isnan(del)
    ierr = 2;
    fprintf('Steepest Descent: Desired tolerance not reached, NaN encountered\n');
else
    ierr = 0;
    fprintf('Steepest descent: Desired tolerance reached in %d iterations\n',cnt);
end
end