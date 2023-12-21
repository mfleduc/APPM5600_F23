clear variables; close all;
I = [0,1];
ns = 2.^(1:9);
fVals = zeros(3,length(ns));
fh = @(x)-4*xlogx(x) ;
for ii = 1:length(ns)
    fVals(1,ii) = compositeIntegration(fh, ns(ii), I,'trapezoid');
    fVals(2,ii) = compositeIntegration(fh, ns(ii), I,'simpson');
    fVals(3,ii) = compositeIntegration(fh, ns(ii), I,'midpoint');
end

figure;
loglog( ns, abs(fVals-1) )
legend( 'Trapezoid', 'Simpsons Rule', 'Midpoint rule' )
grid on
xlabel('n')
ylabel('Error')

figure;
semilogx( ns, abs(fVals-1).*ns.^2 )
legend( 'Trapezoid', 'Simpsons Rule', 'Midpoint rule' )
grid on
xlabel('n')
ylabel('Error*n^2')

