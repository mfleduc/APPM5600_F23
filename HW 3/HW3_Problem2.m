%% HW 3 Problem 2
addpath( '../HW 2')
t = 60*(3600*24);%Time
alpha = 0.138e-6;%Thermal conductivity
Ti = 20;%Initial temperature
Ts = -15;%Snap temperature
tol = 1e-13;%Tolerance, relative error
interval = [0,1];%The interval [0,\bar{x}]
fh = @(x)getTemperature(x, t, alpha, Ti, Ts);
fph = @(x)Tprime(x, t, alpha, Ti, Ts);
maxIters = 100;
%% Run the codes
[root_bisect,cnt, errMsg] = bisection(fh, [0,1], tol, maxIters);
fprintf('Bisection: Error Message %d, Approximate root %.5f\n Next: Newton with x0=0.01\n',errMsg,root_bisect);
[root_n1,errMsgN1,x1] = newtonsMethod(fh, fph, 0.01, tol, maxIters);
fprintf('Newtons Method, initial point 0.01: Error Message %d, Approximate root %.5f\n Next: Newton with x0=1\n',errMsgN1,root_n1);
[root_n2,errMsgN2,x2] = newtonsMethod(fh, fph, interval(2), tol, maxIters);
fprintf('Newtons Method, initial point %d: Error Message %d, Approximate root %.5f\n',interval(2),errMsgN2,root_n2);
%% Speed test
nRuns = 100;
n001Times = zeros(1, nRuns);
n1Times = zeros(1, nRuns);
bisectTimes = zeros(1, nRuns);
for ii = 1:nRuns
    tic;
    [root_bisect,cnt, errMsg] = bisection(fh, [0,1], tol, maxIters);
    bisectTimes(ii)=toc;
    tic;
    [root_n1,errMsgN1,x1] = newtonsMethod(fh, fph, 0.01, tol, maxIters);
    n001Times(ii) = toc;
    tic;
    [root_n2,errMsgN2,x2] = newtonsMethod(fh, fph, interval(2), tol, maxIters);
    n1Times(ii) = toc;
end
