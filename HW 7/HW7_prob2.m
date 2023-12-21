%% HW 7 Problem 2
close all;clear variables;
rng(35683758)
n = 500;
A = 0.5*eye(n);
for ii = 1:n-1
    A(ii,ii+1:end) = 2*(rand([1,n-ii])-0.5);
end

b = randn( [n, 1] );
tol = 10^-10;
A = A+triu(A).';%Building A

tau = [0.01, 0.05, 0.1, 0.2];
residCG = zeros(n,length(tau));
xCG = zeros( n, length(tau) );
evals = zeros(n, length(tau));
xSteep = zeros(size(xCG));
residSteep = zeros(size(residCG));
legendStr = cell(1,length(tau));
itersCG = zeros( n, n,length(tau) );
itersSteep = zeros(n,n,length(tau));
Aj = zeros(n,n,length(tau));
for tt = 1:length(tau)
    legendStr{tt} = sprintf( 'tau = %.2f', tau(tt) );
    mask = (abs(A)<tau(tt)& A~=1);
    Aj(:,:,tt) = A.*mask+eye(n);%Masking A by cutting off large correlations 
    Am = squeeze(Aj(:,:,tt));
    %and making sure not to kill the diagonal
    evals(:,tt) = eig(Am); %CAlculating eigenvalues to determine if SPD: 
    % Will be positive-definite if all eigenvalues are positive
    %% Conjugate gradient
    [xCG(:, tt),errCG, residCG(:,tt), itersCG(:,:,tt)] = ConjGrad( Am, b, tol, n );
    %% Steepest descent
    [xSteep(:,tt),errSteep, residSteep(:,tt), itersSteep(:,:,tt)] = Steepest( Am, b,zeros(size(b)), tol, n);
end

figure;semilogy( 1:n, residCG, 'linewidth', 2) ;
grid on;
legend(legendStr);
xlabel('Iteration number')
ylabel('Relative error')
title('Relative error vs iteration number, conjugate gradient');

figure;semilogy( 1:n, residSteep, 'linewidth', 2) ;
grid on;
legend(legendStr);
xlabel('Iteration number')
ylabel('Relative error')
title('Relative error vs iteration number, Steepest Descent');
xlim([1, n]);

figure;plot(1:n, evals, 'linewidth', 2)
legend(legendStr,'location','southeast');
xlabel('n')
ylabel('\lambda')
title( 'Eigenvalues of Am for each value of tau' )
grid on;

%% Now: Error comparisons
errCGEnergy = zeros(n,length(tau)-1); %Energy norm error
errSteepEnergy = zeros(size(errCGEnergy));
cSteep = zeros(1,3);
cCG = zeros(1,3);
for tt = 1:length(tau)-1
    %Conjugate gradient
    nItersCG = nnz( residCG(:,tt) );
    Am = squeeze(Aj(:,:,tt));
    xstar = Am\b;
    theseItersCG = squeeze(itersCG(:,1:nItersCG,tt));
    ek = theseItersCG-xstar;
    errCGEnergy(1:nItersCG,tt) = diag( (ek')*Am*ek );
    %Steepest descent
    nItersSteep = nnz( residSteep(:,tt) );
    theseItersSteep = squeeze(itersSteep(:,1:nItersSteep,tt));
    ek = theseItersSteep-xstar;
    errSteepEnergy(1:nItersSteep,tt) = diag( (ek')*Am*ek );
    cSteep(tt) = (evals(end,tt)-evals(1,tt))/(evals(end,tt)+evals(1,tt)  );
    k = cond( Am );
    cCG(tt) = (1-sqrt(1/k))/(1+sqrt(1/k));
end
%% Conjugate gradient first
enXstar = xstar'*A*xstar ;
figure;
semilogy( 1:n, enXstar*2*cCG(1).^(1:n), 'r--', 'linewidth', 2 ,'Displayname', 'Error bound, tau=0.01')
grid on;hold on;
semilogy( 1:n, enXstar*2*cCG(2).^(1:n), 'b--', 'linewidth', 2 ,'Displayname', 'Error bound, tau=0.05')
semilogy( 1:n, enXstar*2*cCG(3).^(1:n), 'k--', 'linewidth', 2 ,'Displayname', 'Error bound, tau=0.1')

semilogy( 1:n, errCGEnergy(:,1),'r', 'linewidth', 2,'Displayname', 'Observed error, tau=0.01');
semilogy( 1:n, errCGEnergy(:,2),'b', 'linewidth', 2,'Displayname', 'Observed error, tau=0.05');
semilogy( 1:n, errCGEnergy(:,3),'k', 'linewidth', 2,'Displayname', 'Observed error, tau=0.1');
xlim([1 35])
legend('location', 'southwest')
xlabel('Number of iterations')
ylabel( 'Error' )
title('Error bounds and observed error, conjugate gradient method');
%% STeepest descent
figure; %Note: Initial guess is 0
semilogy( 1:n, enXstar*cSteep(1).^(1:n), 'r--', 'linewidth', 2 ,'Displayname', 'Error bound, tau=0.01')
grid on;hold on;
semilogy( 1:n, enXstar*cSteep(2).^(1:n), 'b--', 'linewidth', 2 ,'Displayname', 'Error bound, tau=0.05')
semilogy( 1:n, enXstar*cSteep(3).^(2:n+1), 'k--', 'linewidth', 2 ,'Displayname', 'Error bound, tau=0.1')
semilogy( 1:n, errSteepEnergy(:,1),'r', 'linewidth', 2,'Displayname', 'Observed error, tau=0.01');
semilogy( 1:n, errSteepEnergy(:,2),'b', 'linewidth', 2,'Displayname', 'Observed error, tau=0.05');
semilogy( 1:n, errSteepEnergy(:,3),'k', 'linewidth', 2,'Displayname', 'Observed error, tau=0.1');
xlim([1 100])
legend('location', 'southeast')
xlabel('Number of iterations')
ylabel( 'Error' )
title('Error bounds and observed error, Steepest descent method');
ylim([10^-30, 10^5])