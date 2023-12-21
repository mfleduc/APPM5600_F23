%% HW 8 Problem 2
clear variables; close all
%Part ii
f = @(x,y)exp(x).*sin(y);%% The function

xi = [0,0,1,1,2,2].';
yi = [0,2,0,2,1,3].';
fi = f(xi,yi);
A = ones(6,6) ;%Design matrix for Ac=f
A(:,2) = xi;
A(:,3) = yi;
A(:,4)=xi.*yi;
A(:,5) = xi.^2;
A(:,6) = yi.^2;
Ainv = inv(A);
c = Ainv*fi;% Solution
%% Part iii
%Surface plot
x = linspace(-1,3,100);
y = x;
[X,Y] = meshgrid(x,y);
F = f(X,Y);%Actual function
% Interpolant
Finterp = c(1)+c(2)*X+c(3)*Y+c(4)*X.*Y +c(5)*X.^2+c(6)*Y.^2;
% plotting
figure;
surf( X,Y,F )%Actual function
xlabel('x');ylabel('y');
zlabel('f(x,y)')
title('Surface plot of f(x,y)')
%Interpolant
figure;
surf( X,Y,Finterp )%Actual function
xlabel('x');ylabel('y');
zlabel('p(x,y)')
title('Surface plot of the interpolant')