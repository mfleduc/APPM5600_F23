%% HW3 Problem 4
clear variables;close all
format long
fh = @(x)(x-1)^5*exp(x); %f(x)
fph = @(x)(x-1)^5*exp(x)+5*(x-1)^4*exp(x);%f'(x)
tol = 10^(-16);%Tolerance
x0=0;%STarting point
maxIters = 1000;%Maximum number of iterations
[root_newton, err_newton, x_newton] = newtonsMethod(fh, fph, x0, tol, maxIters);
figure; semilogy(1:length(x_newton), abs(x_newton-1) )
xlabel('Iteration number');ylabel('Relative error, log10')
grid on
title('Relative error vs iteration number, Newtons method')
saveas(gcf, 'fig_newtons.png');
savefig('fig_newtons.fig')
resN = [[1:10:length(x_newton),length(x_newton)];[x_newton(1:10:end),x_newton(end)];...
    abs([x_newton(1:10:end),x_newton(end)]-1/1) ].';
tableN = array2table(resN, 'variablenames',{'Iteration','x_k', 'relative error'});

[root_mn, err_mn, x_mn] = newtonsMethod(fh, fph, x0, tol, maxIters, 5);
figure; semilogy(1:length(x_mn), abs(x_mn-1) )
xlabel('Iteration number');ylabel('Relative error, log10')
title('Relative error vs iteration number, modified Newtons method')
grid on
saveas(gcf, 'fig_mod.png');
savefig('fig_mod.fig')
resM = [1:length(x_mn);x_mn;abs( x_mn-1 )/1 ].';
tableM = array2table(resM, 'variablenames',{'Iteration','x_k', 'relative error'});