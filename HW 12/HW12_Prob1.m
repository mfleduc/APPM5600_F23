clear variables; close all;
fx = @(x) (1./(1+x.^2));
I = [-5,5];
%% Part B
nTrap = 1291;
nSimp = 108;
resTrap = compositeIntegration(fx, nTrap,I,'trapezoid');
resSimp = compositeIntegration(fx, nSimp,I,'simpson');
[trueVal4,fcnt4] = quad(fx,I(1),I(2),1e-4);
[trueVal6,fcnt6] = quad(fx,I(1),I(2),1e-6);

fprintf('The difference between compositeIntegration and quad for trapeziod rule is %e \n', resTrap-trueVal4);
fprintf('The difference between compositeIntegration and quad for Simpson''s rule is %e \n', resSimp-trueVal4);
fprintf('quad requires %d evaluations for error  less than 1e-4 and %d for error less than 1e-6\n', fcnt4,fcnt6);
