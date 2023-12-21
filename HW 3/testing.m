close all;clear variables;
fh = @(x) 5*x.*(x.^2-1);
fph = @(x)5*(3*x.^2-1 );
gh = @(x) 5*x.^2.*(x-1);
% fpph = @(x) 6*x;
x0 = -2:0.05:2;
% 
figure;plot(x0, fh(x0))
yline(0)
ylim([-10, 10])
hold on;%plot(x0, fph(x0)./fpph(x0))
grid on;

figure;plot(x0, gh(x0))
yline(0)
ylim([-10, 10])
hold on;%plot(x0, fph(x0)./fpph(x0))
grid on;

% for ii = 1:length(x0)
% [root,errMsg, x] = newtonsMethod(fh, fph, x0(ii), 10^-6, 100);
% figure;plot(x)
% if errMsg ~=0
% keyboard
% end
% end
