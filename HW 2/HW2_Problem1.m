%% HW2_Problem1
fh1 = @(x)(x-5)^9;
fh2 = @(x) f1(x);
interval = [4.8, 5.31];
tol = 1e-4;
maxIters = 1000;
[root1, Niters1, errNum1] = bisection( fh1, interval, tol, maxIters );
[root2, Niters2, errNum2] = bisection( fh2, interval, tol, maxIters );

function y = f1(x)
y = zeros(size(x));
for ii = 0:9
    y = y+nchoosek( 9, ii )*(x).^ii*(-5)^( 9-ii );
end
end