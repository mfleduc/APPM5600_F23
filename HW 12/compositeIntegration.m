function res = compositeIntegration(fx,n,I,rule )
%% Function for doing integration via the composite trapezoid rule
% fx: Function handle for the function to be integrated
% I: 2x1 double [a,b] representing the endpoints of the interval to be
% integrated over
% n: Number of subintervals 
% rule: Lowercase string that switches the algorithm between composite trapezoid rule if strcmp(rule,'trapezoid') and 
%   Simpson's rule if strcmp(rule,'simpson')

dx = ( I(2)-I(1) )/n;
interval = linspace( I(1),I(2), n+1 );

switch rule
    case 'trapezoid'
        weights = fx(interval); %R = dx/2*( f(x0)+f(xn)+2*\sum_{k=1}^{n-1} f(xk) )
        weights(2:end-1) = 2*weights(2:end-1);
        res = dx/2*sum(weights);
    case 'simpson'
        weights = fx(interval); %R = dx/2*( f(x0)+f(xn)+2*\sum_{k=1}^{n-1} f(xk) )
        weights(2:end-1) = 2*weights(2:end-1);
        weights(2:2:end-1) = 2*weights(2:2:end-1);
        res = dx/3*sum(weights);
    case 'midpoint'
        midpts = interval(1:end-1)+dx/2;
        weights = fx(midpts);
        res = dx*sum(weights);
    otherwise
        error('Integration rule not defined')
end