function out = InterpHermite(x, fx, fpx, xq)
%% Hermite interpolation
%
xq = reshape( xq, 1, [] );
x = reshape(x,[],1);
fx = reshape(fx,[],1);
fpx = reshape(fpx,[],1);
hVals = zeros(size(xq));
kVals = zeros(size(xq));
for ii = 1:length(x)
    mask = true(size(x));
    mask(ii) = false;
    Li = prod(xq-x(mask))/prod(x(ii)-x(mask));%i-th lagrange polynomial
    ki = (xq-x(ii)).*Li.^2;
    Lip = sum( 1./(x(ii)-x(mask)) );%Derivative of the i-th Lagrange polynomial, evaluated at the interpolation point
    
    hVals = hVals+fx(ii)*(1-2*Lip.*(xq-x(ii))).*Li.^2;
    kVals = kVals+ki*fpx(ii);
end




out = hVals+kVals;
end