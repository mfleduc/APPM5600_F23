%% HW 7 Problem 4
F = @(x,y)[3*x^2+4*y^2-1;y^3-8*x^3-1];%Function

x = [-0.5, 0.25].';

G = [[0.016,-0.17];[0.52,-0.26]];

cnt = 0;
maxIters = 10000;
tol = 10^-7;
del = 1+tol;
while cnt<maxIters && del>tol
   cnt = cnt+1;
   fx = F(x(1),x(2));
   x1 = x-G*fx;
   del = norm(x1-x)/norm(x);
   x = x1;
end