clear variables;close all
f = @(x)(1+x.^2).^(-1);
fp = @(x) -2*x./(1+x.^2).^2;
ns = [5,10,15,20];
finerGrid = linspace(-5,5,1e2+1);h = diff(finerGrid);
twoNorms = zeros( length(ns),8 ); %Record the 2-norms of the differences
supNorms = zeros( length(ns),8 );% Record the sup-norms of the differences
fActual = f( finerGrid );

%% First: Evenly spaced points
for ii = 1:length(ns)
    pts = linspace(-5,5,ns(ii));%The known data points
    fvals = f(pts);%f
    fpvals = fp(pts);%f'
    lagrangeInterp = InterpLagrange(pts, fvals, finerGrid);%Lagrange
    twoNorms(ii,1) = sqrt(h(1))*norm( lagrangeInterp-fActual );%Norms of differences
    supNorms(ii,1)=max(abs(lagrangeInterp-fActual ));
    hermiteInterp = InterpHermite(pts, fvals, fpvals, finerGrid);%Hermite
    twoNorms(ii,2) = sqrt(h(1))*norm( hermiteInterp-fActual );%Norms of differences
    supNorms(ii,2)=max(abs(hermiteInterp-fActual ));
    natSpline = InterpCubicSpline( pts,fvals,finerGrid, 'natural' );%Natural spline
    twoNorms(ii,3) = sqrt(h(1))*norm( natSpline-fActual );%Norms of differences
    supNorms(ii,3)=max(abs(natSpline-fActual ));
    clampedSpline = InterpCubicSpline( pts,fvals,finerGrid, 'clamped',fp([pts(1),pts(end)]) );%Clamped spline
    twoNorms(ii,4) = sqrt(h(1))*norm( clampedSpline-fActual );%Norms of differences
    supNorms(ii,4)=max(abs(clampedSpline-fActual ));
    
    figure;
    plot( finerGrid, f(finerGrid),'k','linewidth', 2 );
    hold on;plot(finerGrid, lagrangeInterp);
    plot(finerGrid, hermiteInterp)
    plot(finerGrid,natSpline)
    plot(finerGrid,clampedSpline)
    plot(pts,fvals,'r.', 'markersize', 20)
    legend('f(x)', 'Lagrange', 'Hermite', 'Natural cubic', 'Clamped cubic','knots')
    ylim([-1.5 1.5])
    xlim([-5.5,5.5])
    grid on;
    xlabel('x')
    ylabel('y')
    title(sprintf('%d equispaced points',ns(ii)));
    saveas(gcf,sprintf('%d equispaced pts.png',ns(ii)))
    savefig(gcf,sprintf('%d equispaced pts.fig',ns(ii)))
end
%% Second: Chebyshev spacing
for ii = 1:length(ns)%First: Evenly spaced points
    pts = -5*cos(pi/2*(2*(1:ns(ii))-1)/ns(ii));%The known data points
    fvals = f(pts);%f
    fpvals = fp(pts);%f'
    lagrangeInterp = InterpLagrange(pts, fvals, finerGrid);%Lagrange
    twoNorms(ii,5) = sqrt(h(1))*norm( lagrangeInterp-fActual );%Norms of differences
    supNorms(ii,5)=max(abs(lagrangeInterp-fActual ));
    hermiteInterp = InterpHermite(pts, fvals, fpvals, finerGrid);%Hermite
    twoNorms(ii,6) = sqrt(h(1))*norm( hermiteInterp-fActual );%Norms of differences
    supNorms(ii,6)=max(abs(hermiteInterp-fActual ));
    natSpline = InterpCubicSpline( pts,fvals,finerGrid, 'natural' );%Natural spline
    twoNorms(ii,7) = sqrt(h(1))*norm( natSpline-fActual );%Norms of differences
    supNorms(ii,7)=max(abs(natSpline-fActual ));
    clampedSpline = InterpCubicSpline( pts,fvals,finerGrid, 'clamped',fp([pts(1),pts(end)]) );%Clamped spline
    twoNorms(ii,8) = sqrt(h(1))*norm( clampedSpline-fActual );%Norms of differences
    supNorms(ii,8)=max(abs(clampedSpline-fActual ));
    
    figure;
    plot( finerGrid, f(finerGrid),'k','linewidth', 2 );
    hold on;plot(finerGrid, lagrangeInterp);
    plot(finerGrid, hermiteInterp)
    plot(finerGrid,natSpline)
    plot(finerGrid,clampedSpline)
    plot(pts,fvals,'r.', 'markersize', 20)
    legend('f(x)', 'Lagrange', 'Hermite', 'Natural cubic', 'Clamped cubic','knots')
    ylim([-1.5 1.5])
    xlim([-5.5,5.5])
    grid on;
    xlabel('x')
    ylabel('y')
    title(sprintf('%d Chebyshev points',ns(ii)));
    saveas(gcf,sprintf('%d Chebyshev pts.png',ns(ii)))
    savefig(gcf,sprintf('%d Chebyshev pts.fig',ns(ii)))
end