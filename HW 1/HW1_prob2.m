%% HW 1 Problem 2
close all; clear variables

x = 1.92:0.001:2.08;
%% Part 1
y1 = f1(x) ;
figure;
plot( x, y1 )
grid on;
xlabel('x')
ylabel('f(x)')
title('Evaluating f(x) expanded')
xlim([1.92, 2.08])
savefig('q2p1.fig')
saveas(gcf, 'q2p1.png')
%% Part 2
f2 = (x-2).^9;
figure;
plot( x, f2 )
grid on
xlabel('x')
ylabel('f(x)')
title('Evaluating f(x) factored')
xlim([1.92, 2.08])
savefig('q2p2.fig')
saveas(gcf, 'q2p2.png')
function y = f1(x)
y = zeros(size(x));
for ii = 0:9
    y = y+nchoosek( 9, ii )*(x).^ii*(-2)^( 9-ii );
end
end