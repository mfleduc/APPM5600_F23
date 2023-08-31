%% Homework 1 Problem 3
clear variables; close all;
%% Part 2
delta = 10.^[-15:0];
x = [0,pi/4, pi/3, pi/2];
figure;
sgtitle( 'Difference between the two functions for x=0,pi/4, pi/3, pi/2' )

for ii = 1:length( x )
    subplot( 2, 2,ii )
    loglog( delta, abs( f1(x(ii),delta )-f2(x(ii), delta )))
    grid on
    xlabel('Delta')
    ylabel('f1(x)-f2(x)')
    title(sprintf('x=%.3f', x(ii)))
end
savefig( 'q3p2.fig' )
saveas(gcf, 'q3p2.png')
%% Part 3
figure;
sgtitle( 'Difference between the two functions for x=0,pi/4, pi/3, pi/2' )

for ii = 1:length( x )
    subplot( 2, 2,ii )
    loglog( delta, abs( f1(x(ii),delta )-f3(x(ii), delta )))
    grid on
    xlabel('Delta')
    ylabel('f1(x)-f3(x)')
    title(sprintf('x=%.3f', x(ii)))
end
savefig( 'q3p3.fig' )
saveas(gcf, 'q3p3.png')


%%%%%%%%%%%%%%%5
function y = f1(x, delta)%%%%Original version of the function
y = cos(x+delta)-cos(x);
end
function y = f2(x, delta)%%%%%New implementation
y = -2*sin(x+delta/2).*sin(delta/2);
end
function y = f3(x, delta)%%%%%%% Taylor expansion
y = -delta*sin(x)+delta.^2.*(-cos((x+delta)/2))/2;
end