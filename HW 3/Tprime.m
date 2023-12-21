function dTdx = Tprime(x, t, alpha, Ti, Ts)


delta = Ti-Ts;
dTdx = 1/sqrt(pi*alpha*t)*exp(-x^2/(4*alpha*t))*delta;
end