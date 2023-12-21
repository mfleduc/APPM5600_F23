function out = InterpCubicSpline( x,fx,xq,bc,varargin )
%% SPline interpolation
%BC: Boundary condition. clamped, natural, or periodic
xq = reshape( xq, 1, [] );
x = reshape(x,[],1);
n = length(x)-1;
h = diff(x);
fx = reshape(fx,[],1);
% fpx = reshape(fpx,[],1);
% fppx = reshape(fppx,[],1);
out = zeros(size(xq));
RHS = zeros(n-1,1);
Vdiag = zeros(n-1,1);
Vband = zeros( n-2,1 );
for ii = 2:n
    RHS(ii-1) = (fx(ii+1)-fx(ii))/h(ii)-(fx(ii)-fx(ii-1))/h(ii-1);%Right hand side of Vm=y
    Vdiag(ii-1) = 1/6*(h(ii)+h(ii-1));
    if ii>2
        Vband(ii-2) = h(ii-1)/6;
    end
end
V = diag(Vdiag);
for ii = 1:n-2
    V(ii,ii+1)=Vband(ii);
end
V = V+V';
Vi = inv(V);
M = Vi*RHS;
switch bc
    case 'clamped'
        fpx = varargin{1};%fpx is a vector with two entries f'(x0) and f'(xn)
        M0 = 3/h(1)*( -fpx(1)+(fx(2)-fx(1))/h(1)-h(1)/6*M(1) );
        Mn = 3/h(end)*( -fpx(end)+(fx(end)-fx(end-1))/h(end)-h(end)/6*M(end) );
        M=[M0;M;Mn];
    case 'natural' %Knots are x0,...,xn so there are n+1 points
        
        M = [0;M;0];    
    case 'periodic'
        error('Periodic boundary conditions are not implemented\n')
    otherwise
        error('Given boundary condition not recognized\n')
end
%% Building the spline
S = zeros( n, length(xq) );
        
for j=1:n
    mask = x(j+1)>xq & x(j)<=xq;
    if j==n
        mask(end)=true;
    end
    
    Cj = fx(j)/h(j)-h(j)*M(j)/6;
    Dj = fx(j+1)/h(j)-h(j)*M(j+1)/6;
    S(j,mask) = 1/(6*h(j))*(M(j)*(x(j+1)-xq(mask)).^3+M(j+1)*(xq(mask)-x(j)).^3)...
        +Cj*(x(j+1)-xq(mask))+Dj*(xq(mask)-x(j));
end
out = sum( S );
end