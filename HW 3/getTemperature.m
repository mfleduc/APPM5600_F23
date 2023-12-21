function temp = getTemperature( x, t, alpha, Ti, Ts )
%% Function for calculating the ground temperature 
%Inputs:
%   x: depth, m, positive down
%   t, time after cold snap, seconds
%   alpha: thermal conductivity, m^2/s
%   Ti: Temperature before cold snap, deg C
%   Ts: Temperature during cold snap, deg C
%Outputs: 
%   temp: Soil temperature at input parameters, deg C
deltaT = Ti-Ts;
temp = erf( x./(2*sqrt( alpha*t ) ))*deltaT+Ts;
end