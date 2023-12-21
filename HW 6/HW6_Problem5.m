%%
clear variables; close all;

A = diag(2*ones([1,6]));%Building A
%Start with a 2 on the diagonal so we can add L+L^T
for ii = 1:5
    A(ii,ii+1)=-1;
    if ii<=3
        A(ii,ii+3)=-1;
    end
end
A = A+A.';
b = [2,1,2,2,1,2].';
x = ones(6,1);

tol = 1e-7;
maxIter = 10000;
%% Gauss-Jacobi method
[xJ,errJ,cJ,normJ] = stationaryIters(A,b,x,tol, maxIter, 'jacobi');
fprintf('Jacobi finished in %d iterations, operator norm is %.4f\n',cJ,normJ);

%% Gauss-Seidel
[xG,errG,cG,normG] = stationaryIters(A,b,x,tol, maxIter, 'gauss-seidel');
fprintf('Gauss-Seidel finished in %d iterations, operator norm is %.4f\n',cG,normG);

%% SOR
w = 1.6735;
[xS,errS,cS,normS] = stationaryIters(A,b,x,tol, maxIter, 'SOR', w);
fprintf('SOR finished in %d iterations, operator norm is %.4f\n\n',cS,normS);
fprintf('log10(Error approximations):\n');
fprintf('Jacobi: %.6f\n', log10(errJ));
fprintf('Gauss-Seidel: %.6f\n', log10(errG));
fprintf('SOR: %.6f\n', log10(errS));

nws = 51;%Number of w values to test
ws = linspace(w-1, w+1, nws);
solns = zeros( 6, nws );
errs = zeros(1,nws);
iters = zeros(1,nws);
for ii = 1:nws
    [solns(:,ii),tmp1,iters(ii),tmp3,errs(ii)] = stationaryIters(A,b,x,tol, maxIter, 'SOR', ws(ii));
end

figure;semilogy(ws(errs==0), iters(errs==0), 'b', 'linewidth',2)
xlabel('w')
ylabel('Number of iterations (Nmax = 10^4)')
title('Number of iterations needed until convergence for different w')
grid on
saveas( gcf, 'sor_wvsiters.png' );
savefig( gcf, 'sor_wvsiters.fig'  )



