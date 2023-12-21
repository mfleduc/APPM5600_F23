function out = InterpLagrange(x,fx, xq)
%% Lagrange interpolating polynomial
% Takes data (x,fx) and interpolates to the query points xq
xq = reshape(xq, 1,[]);
x = reshape(x,[],1);
fx = reshape(fx, [],1);
out = zeros(size(xq));
for ii = 1:length(x)
    mask = true(size(x));
    mask(ii) = false;
    denom = prod( x(ii)-x(mask) );
    num = prod(xq-x(mask));
    out = out+num/denom*fx(ii);
end
end